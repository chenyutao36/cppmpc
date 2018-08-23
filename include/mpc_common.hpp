#ifndef MPC_COMMON_H_
#define MPC_COMMON_H_

typedef struct{
    int nx;        // No. of states
    int nu;        // No. of controls
    int ny;        // No. of cost terms
    int nyN;       // No. of cost terms in the terminal stage
    int np;        // No. of on-line parameters
    int nbx;       // No. of state box constraints
    int nbu;       // No. of control box constraints
    int nbg;       // No. of general constraints
    int nbgN;      // No. of general constraints in the terminal stage
    int N;         // No. of shooting points
    int N2;        // No. of shooting points after partial condensing
    int *nbx_idx;  // index of state box constraints (0-based)
    int *nbu_idx;  // index of control box constraints (0-based)
}model_size;

void regularization(int n, double *A, double reg);

#endif