
#include "mpc_common.hpp"


void regularization(int n, double *A, double reg){
    int i;
    for (i=0;i<n;i++)
        if (A[i*n+i]<reg)
            A[i*n+i] = reg;
}