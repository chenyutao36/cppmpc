#ifndef FULL_CONDENSING_H_
#define FULL_CONDENSING_H_

#include "mpc_common.hpp"
#include "qp_generation.hpp"

typedef struct{
    Matrix<double, Dynamic, Dynamic, RowMajor> Hc;
    Matrix<double, Dynamic, Dynamic, RowMajor> Cc;
    VectorXd gc;
    VectorXd lcc;
    VectorXd ucc;

    MatrixXd G;
    VectorXd L;
    MatrixXd W;
    VectorXd w;
}full_condensing_workspace;

void full_condensing_workspace_init(model_size& size, full_condensing_workspace& cond_work);

void full_condensing(model_size& size, full_condensing_workspace& cond_work,
    qp_in& in, qp_problem& qp, const VectorXd& x0);

#endif