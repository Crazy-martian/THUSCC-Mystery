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