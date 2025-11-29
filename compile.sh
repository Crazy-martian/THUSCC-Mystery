source ~/.bashrc
conda activate mystery_env

export CFITSIO_INSTALL_DIR=$PWD/cfitsio_local

echo "compiling"

g++ ratio_calculate.cpp -o ratio_calc \
    -std=c++17 \
    # -fopenmp \
    -O3 \
    -I$CFITSIO_INSTALL_DIR/include \    # 告诉 g++ 在哪里查找 fitsio.h
    -L$CFITSIO_INSTALL_DIR/lib \       # 告诉 g++ 在哪里查找 libcfitsio.so
    -lcfitsio \                          # 告诉 g++ 链接 cfitsio 库
    -lm                                # 链接数学库

