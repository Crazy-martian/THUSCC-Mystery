import numpy as np
from astropy.io import fits
import math
from numba import njit, prange
import time
DATA_SIZE_X = 512
DATA_SIZE_Y = 512
CENTER_X = 256.0
CENTER_Y = 256.0
# 网格大小设定为 16 (略大于 PSF_radius 7.5，保证只需查相邻格子)
CELL_SIZE = 16 
GRID_W = DATA_SIZE_X // CELL_SIZE
GRID_H = DATA_SIZE_Y // CELL_SIZE
MAX_EVENTS_PER_CELL = 100  # 假设每个格子最多这么多光子，防止溢出，可视情况调整

@njit(fastmath=True, cache=True)
def precompute_and_build_grid(data):
    # --- Geometry Map (Float32) ---
    # 这一步其实可以用 global memory 缓存，避免每次都算，但这里为了纯净算在耗时里
    # 简化：因为 cos/sin 对称，其实可以优化，但 Numba 这里不是瓶颈
    sigma_major_map = np.empty((DATA_SIZE_X, DATA_SIZE_Y), dtype=np.float32)
    sigma_minor_map = np.empty((DATA_SIZE_X, DATA_SIZE_Y), dtype=np.float32)
    sin_map = np.empty((DATA_SIZE_X, DATA_SIZE_Y), dtype=np.float32)
    cos_map = np.empty((DATA_SIZE_X, DATA_SIZE_Y), dtype=np.float32)
    norm_map = np.empty((DATA_SIZE_X, DATA_SIZE_Y), dtype=np.float32)

    max_dist_inv = 1.0 / 362.039

    for x in range(DATA_SIZE_X):
        dx = CENTER_X - x - 0.5
        for y in range(DATA_SIZE_Y):
            dy = CENTER_Y - y - 0.5
            dist = math.sqrt(dx*dx + dy*dy)
            
            ratio = dist * max_dist_inv
            sigma_minor = 0.5 + ratio * 2.5
            eccentricity = 0.9 * ratio
            
            ecc_sq = eccentricity * eccentricity
            denom = 0.001 if ecc_sq >= 1.0 else math.sqrt(1.0 - ecc_sq)
            sigma_major = sigma_minor / denom

            angle = math.atan2(dy, dx)
            
            sigma_major_map[x, y] = sigma_major
            sigma_minor_map[x, y] = sigma_minor
            sin_map[x, y] = math.sin(angle)
            cos_map[x, y] = math.cos(angle)
            norm_map[x, y] = 1.0 / (6.283185 * sigma_major * sigma_minor)

    # --- Grid Indexing ---
    n_events = len(data)
    grid_counts = np.zeros((GRID_W, GRID_H), dtype=np.int32)
    grid_indices = np.zeros((GRID_W, GRID_H, MAX_EVENTS_PER_CELL), dtype=np.int32)

    for i in range(n_events):
        ev_x = data[i, 0]
        ev_y = data[i, 1]
        gx = int(ev_x) // CELL_SIZE
        gy = int(ev_y) // CELL_SIZE
        
        if 0 <= gx < GRID_W and 0 <= gy < GRID_H:
            c = grid_counts[gx, gy]
            if c < MAX_EVENTS_PER_CELL:
                grid_indices[gx, gy, c] = i
                grid_counts[gx, gy] += 1

    return sigma_major_map, sigma_minor_map, sin_map, cos_map, norm_map, grid_counts, grid_indices

@njit(parallel=True, fastmath=True, cache=True)
def core_calculation(data, bkg_rate=1e-5, t_exp=1000.0):
    # 1. 预计算和构建网格
    s_maj, s_min, sa, ca, nc, g_cnt, g_idx = precompute_and_build_grid(data)
    
    ratio_grid = np.zeros((DATA_SIZE_X, DATA_SIZE_Y), dtype=np.float32)
    R_grid = np.zeros((DATA_SIZE_X, DATA_SIZE_Y), dtype=np.float32)

    psf_radius = 7.0
    iteration = 5
    bkg_rate_f = np.float32(bkg_rate)
    t_f = np.float32(t_exp)
    inv_bkg = 1.0 / bkg_rate_f

    # 2. 遍历 Grid (而不是遍历所有 Pixel)
    # prange 并行化 grid 遍历
    for gx in prange(GRID_W):
        for gy in range(GRID_H):
            # 检查当前 Grid 以及周围 3x3 邻域是否有光子
            # 如果邻域全是空的，这个 Grid 里的所有像素都不需要计算
            has_events = False
            for ix in range(max(0, gx-1), min(GRID_W, gx+2)):
                for iy in range(max(0, gy-1), min(GRID_H, gy+2)):
                    if g_cnt[ix, iy] > 0:
                        has_events = True
                        break
                if has_events: break
            
            if not has_events:
                continue

            # 3. 如果是活跃 Grid，遍历其中的像素
            x_start = gx * CELL_SIZE
            y_start = gy * CELL_SIZE
            # 边界保护
            x_end = min(x_start + CELL_SIZE, DATA_SIZE_X)
            y_end = min(y_start + CELL_SIZE, DATA_SIZE_Y)

            for x in range(x_start, x_end):
                for y in range(y_start, y_end):
                    
                    # --- Pixel Kernel Start ---
                    # 收集光子 (Gather)
                    # 只需搜索 3x3 的邻居格子
                    local_count = 0
                    local_s = np.zeros(256, dtype=np.float32) # Stack array

                    # Pixel geometric params
                    px_smaj = s_maj[x, y]
                    px_smin = s_min[x, y]
                    px_sa = sa[x, y]
                    px_ca = ca[x, y]
                    px_nc = nc[x, y]
                    
                    # 确定要搜索的 Grid 范围 (当前像素所属 Grid 的邻居)
                    # 注意：这里 x, y 属于 gx, gy，所以只需搜 gx, gy 的邻居
                    for ix in range(max(0, gx-1), min(GRID_W, gx+2)):
                        for iy in range(max(0, gy-1), min(GRID_H, gy+2)):
                            cnt = g_cnt[ix, iy]
                            if cnt == 0: continue
                            
                            for k in range(cnt):
                                idx = g_idx[ix, iy, k]
                                dx_ev = data[idx, 0] - x
                                dy_ev = data[idx, 1] - y
                                
                                # Bounding Box Check
                                if abs(dx_ev) <= psf_radius and abs(dy_ev) <= psf_radius:
                                    # PSF Calc
                                    t1 = (dx_ev * px_ca + dy_ev * px_sa) / px_smaj
                                    t2 = (-dx_ev * px_sa + dy_ev * px_ca) / px_smin
                                    expo = (t1*t1 + t2*t2) * 0.5
                                    val = math.exp(-expo) * px_nc
                                    
                                    if local_count < 256:
                                        local_s[local_count] = val
                                        local_count += 1
                    
                    if local_count == 0:
                        continue

                    # R Iteration
                    R = np.float32(0.05)
                    for _ in range(iteration):
                        temp_R = np.float32(0.0)
                        for k in range(local_count):
                            s = local_s[k]
                            temp_R += (R * s) / ((R * s + bkg_rate_f) * t_f)
                        R = temp_R
                    
                    R_grid[x, y] = R

                    # Log Likelihood
                    log_sum = np.float32(0.0)
                    for k in range(local_count):
                        # term = (R*s + bkg)/bkg = 1 + R*s*inv_bkg
                        log_sum += np.log(1.0 + R * local_s[k] * inv_bkg)
                    
                    ratio_grid[x, y] = log_sum - (t_f * R)
                    # --- Pixel Kernel End ---

    return ratio_grid, R_grid

@njit(fastmath=True, cache=True)
def find_local_maxima(ratio_grid, R_grid, threshold=8, window_size=21, border=10):
    """
    替代 scipy.ndimage.maximum_filter。
    只在大于阈值的点进行邻域检查，性能比 scipy 快得多。
    """
    w_radius = window_size // 2
    sources = []
    
    # 直接遍历，但加上步长或快速跳过小于阈值的区域
    # 这里为了简单直接全遍历，因为 Numba 循环极快
    rows, cols = ratio_grid.shape
    
    for r in range(border, rows - border):
        for c in range(border, cols - border):
            val = ratio_grid[r, c]
            
            if val < threshold:
                continue
            
            # 如果大于阈值，检查是否是局部最大
            is_max = True
            # 手动遍历窗口
            r_min = max(0, r - w_radius)
            r_max = min(rows, r + w_radius + 1)
            c_min = max(0, c - w_radius)
            c_max = min(cols, c + w_radius + 1)
            
            for i in range(r_min, r_max):
                for j in range(c_min, c_max):
                    if i == r and j == c:
                        continue
                    if ratio_grid[i, j] > val: # 严格大于，或者 >= 看来题意
                        is_max = False
                        break
                if not is_max:
                    break
            
            if is_max:
                # 记录结果: x, y, rate. 注意 data index 和坐标的对应
                # 原代码是 data[n, 0] 是 x, 对应 array 的 index 可能是 [x, y] 或 [y, x]
                # 假设 ratio_grid[x, y] 对应物理坐标 x, y
                sources.append((float(r), float(c), R_grid[r, c]))
                
    return sources

if __name__ == '__main__':
    t0 = time.time()

    # 读取数据
    hdu = fits.open('mock_data.fits', memmap=True)
    x_coord = hdu[1].data['x']
    y_coord = hdu[1].data['y']
    data = np.column_stack((x_coord, y_coord)).astype(np.float64)
    data = np.ascontiguousarray(data)

    t1 = time.time()
    print(f"Data IO time: {t1 - t0:.4f}s")

    # 计算
    print("Start calculating map ...")
    ratio_grid, R_grid = core_calculation(data)

    t2 = time.time()
    print(f"Map calculation time: {t2 - t1:.4f}s")

    # 检测
    source_coord_list = find_local_maxima(ratio_grid, R_grid)
    print(f"Detected sources: {len(source_coord_list)}")

    t3 = time.time()
    print(f"Detection time: {t3 - t2:.4f}s")

    # 保存结果
    source_coord_array = np.array(source_coord_list)
    # col1 = fits.Column(name='x', format='I', array=source_coord_array[:, 0])
    # col2 = fits.Column(name='y', format='I', array=source_coord_array[:, 1])
    # col3 = fits.Column(name='countrate', format='D', array=source_coord_array[:, 2])
    # cols = fits.ColDefs([col1, col2, col3])
    # tbhdu = fits.BinTableHDU.from_columns(cols)
    # prihdr = fits.Header()
    # prihdr['COMMENT'] = "This file storages the info of detected sources"
    # prihdu = fits.PrimaryHDU(header=prihdr)
    # hdulist = fits.HDUList([prihdu, tbhdu])
    # hdulist.writeto('detection_info.fits', overwrite=True)

    with fits.open('detection_info.fits', mode='update') as hdul:
        tb_hdu = hdul[1]
        tb_hdu.data['x'] = source_coord_array[:, 0]
        tb_hdu.data['y'] = source_coord_array[:, 1]
        tb_hdu.data['countrate'] = source_coord_array[:, 2]

    t4 = time.time()
    print(f"Data IO time: {t4 - t3:.4f}s")
    print(f"Total time: {t4 - t0:.4f}s")
