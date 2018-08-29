# cppmpc

Installation Guide:

1. Download qpoases-3.2.0 at https://www.coin-or.org/download/source/qpOASES/ and install it linking with lapack and blas. (do not use qpoases-3.2.1 since there are some missing symbols while linking)

2. In make.mk, change QPOASES_PATH to the directory of qpoases-3.2.0

3. Run make
