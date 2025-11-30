  source ~/.bashrc
conda activate mystery_env

export CFITSIO_INSTALL_DIR=$PWD/cfitsio_local

echo "compiling"

g++ ratio_calculate.cpp -o ratio_calc -std=c++17 -fopenmp -O3 -I$CFITSIO_INSTALL_DIR/include -L$CFITSIO_INSTALL_DIR/lib  -lcfitsio -lm  -ffast-math                             

