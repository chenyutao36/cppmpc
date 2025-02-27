# cppmpc

## Installation Guide:

1. Download Eigen and extract (No install needed)

2. Download qpoases-3.2.0 at https://www.coin-or.org/download/source/qpOASES/ (do not use qpoases-3.2.1 since there are some missing symbols while linking)

3. Install qpoases using lapack and blas (in Makefile.linux set REPLACE_LINALG = 0, and set the paths for LIB_BLAS and LIB_LAPACK

4. In make.mk, change QPOASES_PATH to the directory of qpoases-3.2.0, and set the Eigen path.

5. Run make

### To run your own examples, use MATMPC to generate casadi_src.cpp and put it in the /model folder. Create your own main function in the /examples folder and add it in the makefile list.
