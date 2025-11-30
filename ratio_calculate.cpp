#include <iostream>
#include <vector>
#include <string>
#include <memory>    // for std::unique_ptr
#include <stdexcept> // for std::runtime_error
#include <cmath>
#include <chrono>
#include <omp.h>
#include "fitsio.h"

constexpr int DATA_SIZE_X = 512;
constexpr int DATA_SIZE_Y = 512;
constexpr float CENTER_X = 256.0f;
constexpr float CENTER_Y = 256.0f;
constexpr int CELL_SIZE = 16; 
constexpr int GRID_W = DATA_SIZE_X / CELL_SIZE;
constexpr int GRID_H = DATA_SIZE_Y / CELL_SIZE;

constexpr int MAX_EVENTS_PER_CELL = 100; 
constexpr float PSF_RADIUS = 7.0f;       
constexpr int ITERATION = 5;             
constexpr float BKG_RATE = 1e-5f;        
constexpr float T_EXP = 1000.0f;         
constexpr float INV_BKG = 1.0f / BKG_RATE;
constexpr float INITIAL_R = 0.05f;

constexpr float THRESHOLD = 8.0f;      
constexpr int WINDOW_SIZE = 21;        
constexpr int W_RADIUS = WINDOW_SIZE / 2; 
constexpr int BORDER = 10;               

// Helper function to get 2D index from flattened 1D array
inline int idx(int x, int y) {
    return x * DATA_SIZE_Y + y;
}

struct Event {
    double x;
    double y;
};

struct Source {
    float x;
    float y;
    float rate;
};

void check_fits_status(int status) {
    if (status) {
        char err_text[FLEN_ERRMSG];
        fits_read_errmsg(err_text);
        throw std::runtime_error(std::string("CFITSIO Error: ") + err_text);
    }
}



/**
Performs the core likelihood calculation. This function ports both precompute_and_build_grid and core_calculation from the Python script.
data       Vector of all photon events.
ratio_grid Output pointer for the likelihood map (must be pre-allocated).
R_grid     Output pointer for the rate map (must be pre-allocated).
 */
void core_calculation(const std::vector<Event>& data, float* ratio_grid, float* R_grid) {
    
    // --- 1. Precomputation (from precompute_and_build_grid) ---
    
    // --- 1a. Geometry Maps ---
    // Allocate memory for PSF geometry maps
    auto s_maj_map_ptr = std::make_unique<float[]>(DATA_SIZE_X * DATA_SIZE_Y);
    auto s_min_map_ptr = std::make_unique<float[]>(DATA_SIZE_X * DATA_SIZE_Y);
    auto sin_map_ptr   = std::make_unique<float[]>(DATA_SIZE_X * DATA_SIZE_Y);
    auto cos_map_ptr   = std::make_unique<float[]>(DATA_SIZE_X * DATA_SIZE_Y);
    auto norm_map_ptr  = std::make_unique<float[]>(DATA_SIZE_X * DATA_SIZE_Y);
    
    // Get raw pointers for easier access
    float* s_maj = s_maj_map_ptr.get();
    float* s_min = s_min_map_ptr.get();
    float* sa    = sin_map_ptr.get();
    float* ca    = cos_map_ptr.get();
    float* nc    = norm_map_ptr.get();

    const float max_dist_inv = 1.0f / 362.039f; 
    const float pi = 3.1415926535f;

    // This loop is not parallelized in Python
    for (int x = 0; x < DATA_SIZE_X; ++x) {
        float dx = CENTER_X - x - 0.5f; 
        for (int y = 0; y < DATA_SIZE_Y; ++y) {
            float dy = CENTER_Y - y - 0.5f; 
            float dist = std::sqrt(dx * dx + dy * dy);
            
            float ratio = dist * max_dist_inv;
            float sigma_minor = 0.5f + ratio * 2.5f; 
            float eccentricity = 0.9f * ratio;     
            
            float ecc_sq = eccentricity * eccentricity;
            float denom = (ecc_sq >= 1.0f) ? 0.001f : std::sqrt(1.0f - ecc_sq); 
            float sigma_major = sigma_minor / denom;

            float angle = std::atan2(dy, dx); 
            
            int index = idx(x, y);
            s_maj[index] = sigma_major;
            s_min[index] = sigma_minor;
            sa[index]    = std::sin(angle);
            ca[index]    = std::cos(angle);
            nc[index]    = 1.0f / (2.0f * pi * sigma_major * sigma_minor); 
        }
    }

    // --- 1b. Grid Indexing ---
    auto grid_counts_ptr = std::make_unique<int[]>(GRID_W * GRID_H);
    // This is a 3D flattened array: grid[gx][gy][event_index]
    auto grid_indices_ptr = std::make_unique<int[]>(GRID_W * GRID_H * MAX_EVENTS_PER_CELL);
    
    int* g_cnt = grid_counts_ptr.get();
    int* g_idx = grid_indices_ptr.get();

    // Initialize grid counts to zero
    std::fill_n(g_cnt, GRID_W * GRID_H, 0);

    long n_events = data.size();
    for (int i = 0; i < n_events; ++i) {
        int gx = static_cast<int>(data[i].x) / CELL_SIZE; 
        int gy = static_cast<int>(data[i].y) / CELL_SIZE;
        
        if (gx >= 0 && gx < GRID_W && gy >= 0 && gy < GRID_H) {
            int grid_index = gx * GRID_H + gy; // 2D index for grid_counts
            int c = g_cnt[grid_index];
            if (c < MAX_EVENTS_PER_CELL) {
                // 3D index for grid_indices
                g_idx[grid_index * MAX_EVENTS_PER_CELL + c] = i;
                g_cnt[grid_index]++; 
            }
        }
    }
    
    // --- 2. Core Calculation (Parallelized) ---
    
    // Initialize output grids to zero
    std::fill_n(ratio_grid, DATA_SIZE_X * DATA_SIZE_Y, 0.0f);
    std::fill_n(R_grid, DATA_SIZE_X * DATA_SIZE_Y, 0.0f);

    // This is the C++/OpenMP equivalent of 'prange(GRID_W)'
    #pragma omp parallel for
    for (int gx = 0; gx < GRID_W; ++gx) {
        
        // This array is local to each thread, simulating the 'stack array'
        float local_s[256]; 

        for (int gy = 0; gy < GRID_H; ++gy) {
            
            // --- Active Grid Check ---
            //
            bool has_events = false;
            for (int ix = std::max(0, gx - 1); ix < std::min(GRID_W, gx + 2); ++ix) {
                for (int iy = std::max(0, gy - 1); iy < std::min(GRID_H, gy + 2); ++iy) {
                    if (g_cnt[ix * GRID_H + iy] > 0) {
                        has_events = true;
                        break;
                    }
                }
                if (has_events) break;
            }
            if (!has_events) continue; 

            // --- Pixel Loop ---
            //
            int x_start = gx * CELL_SIZE;
            int y_start = gy * CELL_SIZE;
            int x_end = std::min(x_start + CELL_SIZE, DATA_SIZE_X);
            int y_end = std::min(y_start + CELL_SIZE, DATA_SIZE_Y);

            for (int x = x_start; x < x_end; ++x) {
                for (int y = y_start; y < y_end; ++y) {
                    
                    int local_count = 0;
                    int p_idx = idx(x, y); // Pixel index
                    
                    // --- Gather Photons from 3x3 Grid ---
                    for (int ix = std::max(0, gx - 1); ix < std::min(GRID_W, gx + 2); ++ix) {
                        for (int iy = std::max(0, gy - 1); iy < std::min(GRID_H, gy + 2); ++iy) {
                            
                            int grid_index = ix * GRID_H + iy;
                            int cnt = g_cnt[grid_index];
                            if (cnt == 0) continue; 

                            for (int k = 0; k < cnt; ++k) {
                                int event_idx = g_idx[grid_index * MAX_EVENTS_PER_CELL + k];
                                float dx_ev = data[event_idx].x - x;
                                float dy_ev = data[event_idx].y - y;
                                
                                // Bounding Box Check
                                if (std::abs(dx_ev) <= PSF_RADIUS && std::abs(dy_ev) <= PSF_RADIUS) {
                                    // PSF Calc
                                    float t1 = (dx_ev * ca[p_idx] + dy_ev * sa[p_idx]) / s_maj[p_idx];
                                    float t2 = (-dx_ev * sa[p_idx] + dy_ev * ca[p_idx]) / s_min[p_idx];
                                    float expo = (t1 * t1 + t2 * t2) * 0.5f;
                                    float val = std::exp(-expo) * nc[p_idx];
                                    
                                    if (local_count < 256) {
                                        local_s[local_count] = val; 
                                        local_count++;
                                    }
                                }
                            }
                        }
                    } // end 3x3 grid search
                    
                    if (local_count == 0) continue; 

                    // --- R Iteration ---
                    float R = INITIAL_R;
                    for (int iter = 0; iter < ITERATION; ++iter) {
                        float temp_R = 0.0f;
                        for (int k = 0; k < local_count; ++k) {
                            float s = local_s[k];
                            temp_R += (R * s) / ((R * s + BKG_RATE) * T_EXP);
                        }
                        R = temp_R;
                    }
                    R_grid[p_idx] = R; 

                    // --- Log Likelihood ---
                    float log_sum = 0.0f;
                    for (int k = 0; k < local_count; ++k) {
                        log_sum += std::log(1.0f + R * local_s[k] * INV_BKG);
                    }
                    ratio_grid[p_idx] = log_sum - (T_EXP * R);

                } // end y
            } // end x
        } // end gy
    } // end gx (parallel loop)
}

std::vector<Source> find_local_maxima(const float* ratio_grid, const float* R_grid) {
    std::vector<Source> sources;

    // This loop can also be parallelized, as each pixel is independent
    // We add 'schedule(dynamic)' because sources are sparse, so some
    // threads might finish their rows much faster than others.
    #pragma omp parallel for schedule(dynamic)
    for (int r = BORDER; r < DATA_SIZE_X - BORDER; ++r) {
        for (int c = BORDER; c < DATA_SIZE_Y - BORDER; ++c) {
            
            int p_idx = idx(r, c); // Use our helper: r * DATA_SIZE_Y + c
            float val = ratio_grid[p_idx];

            // 1. Check threshold
            if (val < THRESHOLD) {
                continue;
            }

            // 2. Check if it's a local maximum
            bool is_max = true;
            for (int i = r - W_RADIUS; i <= r + W_RADIUS; ++i) {
                for (int j = c - W_RADIUS; j <= c + W_RADIUS; ++j) {
                    if (i == r && j == c) {
                        continue;
                    }
                    // Check if neighbor is greater
                    if (ratio_grid[idx(i, j)] > val) { 
                        is_max = false;
                        break;
                    }
                }
                if (!is_max) {
                    break;
                }
            }

            // 3. If it is, add to list
            if (is_max) {
                // We must use 'omp critical' to safely add to the
                // shared 'sources' vector from multiple threads.
                #pragma omp critical
                {
                    sources.push_back({static_cast<float>(r), 
                                       static_cast<float>(c), 
                                       R_grid[p_idx]});
                }
            }
        } // end for c
    } // end for r (parallel)

    return sources;
}

void write_fits_results(const std::string& filename, const std::vector<Source>& sources) {
    fitsfile* fptr = nullptr; // C-style file pointer
    int status = 0;   // CFITSIO status variable
    
    long n_sources = sources.size();
    if (n_sources == 0) {
        std::cout << "No sources to write, skipping FITS update." << std::endl;
        return;
    }

    try {
        // 1. Open the FITS file in READWRITE mode
        fits_open_file(&fptr, filename.c_str(), READWRITE, &status);
        check_fits_status(status);

        // 2. Move to the 2nd HDU (binary table)
        int hdutype;
        fits_movabs_hdu(fptr, 2, &hdutype, &status);
        check_fits_status(status);

        // 3. Get column numbers
        int x_col, y_col, rate_col;
        fits_get_colnum(fptr, CASEINSEN, (char*)"x", &x_col, &status);
        check_fits_status(status);
        fits_get_colnum(fptr, CASEINSEN, (char*)"y", &y_col, &status);
        check_fits_status(status);
        fits_get_colnum(fptr, CASEINSEN, (char*)"countrate", &rate_col, &status);
        check_fits_status(status);

        // 4. Prepare C-style arrays from the std::vector<Source>
        // must match the data types in the FITS file.
        // 'x' and 'y' are format 'I' (long integer in C)
        // 'countrate' is format 'D' (double in C)
        std::unique_ptr<long[]> x_data(new long[n_sources]);
        std::unique_ptr<long[]> y_data(new long[n_sources]);
        std::unique_ptr<double[]> rate_data(new double[n_sources]);

        for (long i = 0; i < n_sources; ++i) {
            // Cast our float results to the required FITS types
            x_data[i]    = static_cast<long>(sources[i].x);
            y_data[i]    = static_cast<long>(sources[i].y);
            rate_data[i] = static_cast<double>(sources[i].rate);
        }

        // 5. Write the data to the columns
        long firstrow = 1;
        long firstelem = 1;
        fits_write_col(fptr, TLONG,   x_col,    firstrow, firstelem, n_sources, x_data.get(),    &status);
        check_fits_status(status);
        fits_write_col(fptr, TLONG,   y_col,    firstrow, firstelem, n_sources, y_data.get(),    &status);
        check_fits_status(status);
        fits_write_col(fptr, TDOUBLE, rate_col, firstrow, firstelem, n_sources, rate_data.get(), &status);
        check_fits_status(status);

        // 6. Close the FITS file
        fits_close_file(fptr, &status);
        check_fits_status(status);
        
        std::cout << "Successfully wrote " << n_sources << " sources to " << filename << std::endl;

    } catch (const std::exception& e) {
        if (fptr) fits_close_file(fptr, &status);
        std::cerr << "Error writing to FITS file: " << e.what() << std::endl;
    }
}

std::vector<Event> read_fits_data(const std::string& filename) {
    fitsfile* fptr; 
    int status = 0;   
    
    try {
        
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);
        check_fits_status(status);
       
        int hdutype;
        fits_movabs_hdu(fptr, 2, &hdutype, &status);
        check_fits_status(status);

        if (hdutype != BINARY_TBL) {
            throw std::runtime_error("FITS file does not contain a binary table at HDU 2.");
        }

        long nrows;
        fits_get_num_rows(fptr, &nrows, &status);
        check_fits_status(status);

       
        int x_col, y_col;
        fits_get_colnum(fptr, CASEINSEN, (char*)"x", &x_col, &status);
        check_fits_status(status);
        fits_get_colnum(fptr, CASEINSEN, (char*)"y", &y_col, &status);
        check_fits_status(status);
     
        std::unique_ptr<double[]> x_data(new double[nrows]);
        std::unique_ptr<double[]> y_data(new double[nrows]);
     
        long firstrow = 1;
        long firstelem = 1;
        double nulval = 0.0; 
        int anynul;      

        fits_read_col(fptr, TDOUBLE, x_col, firstrow, firstelem, nrows, &nulval, x_data.get(), &anynul, &status);
        check_fits_status(status);
        fits_read_col(fptr, TDOUBLE, y_col, firstrow, firstelem, nrows, &nulval, y_data.get(), &anynul, &status);
        check_fits_status(status);

        std::vector<Event> data;
        data.resize(nrows);
        for (long i = 0; i < nrows; ++i) {
            data[i].x = x_data[i];
            data[i].y = y_data[i];
        }

        fits_close_file(fptr, &status);
        check_fits_status(status);
        
        return data;

    } catch (const std::exception& e) {
        if (fptr) fits_close_file(fptr, &status);
        std::cerr << e.what() << std::endl;
        return std::vector<Event>();
    }
}

int main() {
    try {
        auto t0 = std::chrono::high_resolution_clock::now();

        // 1. Read Data
        std::cout << "Reading data..." << std::endl;
        std::vector<Event> data = read_fits_data("mock_data.fits");
        
        if (data.empty()) {
            std::cerr << "Failed to read data. Exiting." << std::endl;
            return 1;
        }
        std::cout << "Successfully read " << data.size() << " events." << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        std::cout << "Data IO time: " 
                  << std::chrono::duration<double>(t1 - t0).count() << "s" << std::endl;

        // 2. Allocate memory for grids
        auto ratio_grid_ptr = std::make_unique<float[]>(DATA_SIZE_X * DATA_SIZE_Y);
        auto R_grid_ptr     = std::make_unique<float[]>(DATA_SIZE_X * DATA_SIZE_Y);
        float* ratio_grid = ratio_grid_ptr.get();
        float* R_grid     = R_grid_ptr.get();

        // 3. Calculate
        std::cout << "Start calculating map..." << std::endl;
        core_calculation(data, ratio_grid, R_grid); 

        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "Map calculation time: " 
                  << std::chrono::duration<double>(t2 - t1).count() << "s" << std::endl;

        // 4. Detect Sources (Stub)
        std::vector<Source> sources = find_local_maxima(ratio_grid, R_grid); 
        std::cout << "Detected sources (stub): " << sources.size() << std::endl;

        auto t3 = std::chrono::high_resolution_clock::now();
        std::cout << "Detection time (stub): " 
                  << std::chrono::duration<double>(t3 - t2).count() << "s" << std::endl;

        // 5. Write Results (Stub)
        write_fits_results("detection_info.fits", sources); 

        auto t4 = std::chrono::high_resolution_clock::now();
        std::cout << "Data IO time (stub): " 
                  << std::chrono::duration<double>(t4 - t3).count() << "s" << std::endl;
        
        std::cout << "Total time: " 
                  << std::chrono::duration<double>(t4 - t0).count() << "s" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "An unexpected error occurred: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}




int __previous_main_for_testing_io() {
    try {
        auto t0 = std::chrono::high_resolution_clock::now();

        std::cout << "Attempting to read 'mock_data.fits'..." << std::endl;
        std::vector<Event> data = read_fits_data("mock_data.fits");

        auto t1 = std::chrono::high_resolution_clock::now();

        if (data.empty()) {
            std::cerr << "Failed to read data." << std::endl;
            return 1;
        }

        std::cout << "Successfully read " << data.size() << " events." << std::endl;
        std::cout << "Data IO time: " << std::chrono::duration<double>(t1 - t0).count() << "s" << std::endl;


        std::cout << "First 5 events:" << std::endl;
        for (int i = 0; i < 5 && i < data.size(); ++i) {
            std::cout << "Event " << i << ": (x=" << data[i].x << ", y=" << data[i].y << ")" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "An unexpected error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

