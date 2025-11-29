import numpy as np
from astropy.io import fits

def score_calc():
    '''
    This function calculates the score of the detection results
    '''
    hdu0 = fits.open('source_info.fits')
    x_coord = hdu0[1].data['x']
    y_coord = hdu0[1].data['y']
    rate = hdu0[1].data['countrate']
    source_data = np.column_stack((x_coord,y_coord,rate))
    hdu0.close()

    hdu1 = fits.open('detection_info.fits')
    x_coord = hdu1[1].data['x']
    y_coord = hdu1[1].data['y']
    rate = hdu1[1].data['countrate']
    det_data = np.column_stack((x_coord,y_coord,rate))
    hdu1.close()

    max_min_distance = 0
    max_delta_rate = 0

    for i in range(12):
        print(f"source rate of {i} is {source_data[i][2]}")
    for i in range(12):
        print(f"det rate of {i} is {det_data[i][2]}")
    for i in range (12):
        det_x = det_data[i][0]
        det_y = det_data[i][1]
        det_rate = det_data[i][2]
        distance = 500
        for j in range(12):
            source_x = source_data[j][0]
            source_y = source_data[j][1]
            source_rate = source_data[j][2]
            distance1 = ((source_x - det_x)**2+(source_y - det_y)**2)**0.5
            if distance1 < distance:
                distance = distance1
                delta_rate = abs(source_rate - det_rate)

        print(f"source #{i}: distance = {distance}, delta rate = {delta_rate}")
        max_min_distance = max(max_min_distance, distance)
        max_delta_rate = max(max_delta_rate, delta_rate)

    print(f"max distance: {max_min_distance}")
    print(f"max delta rate: {max_delta_rate}")

if __name__ == '__main__':
    score_calc()