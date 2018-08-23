#ifndef QP_GENERATION_H_
#define QP_GENERATION_H_

#include <Eigen/Dense>
#include "mpc_common.hpp"

using namespace Eigen;

typedef struct{

    MatrixXd x;
    MatrixXd u; 
    MatrixXd y;
    VectorXd yN;
    MatrixXd W;
    MatrixXd WN;
    MatrixXd p;
    VectorXd lbx;
    VectorXd ubx;
    VectorXd lbu;
    VectorXd ubu;
    VectorXd lbg;
    VectorXd ubg;
    VectorXd lbgN;
    VectorXd ubgN;
    double reg;
}qp_in;

typedef struct{
    MatrixXd Q;
    MatrixXd S;
    MatrixXd R;
    MatrixXd A;
    MatrixXd B;
    MatrixXd a;
    MatrixXd gx;
    MatrixXd gu;
    MatrixXd Cx;
    MatrixXd Cgx;
    MatrixXd CgN;
    MatrixXd Cgu;
    VectorXd lb_u;
    VectorXd ub_u;
    VectorXd lb_x;
    VectorXd ub_x;
    VectorXd lb_g;
    VectorXd ub_g;
}qp_problem;

typedef struct{
    MatrixXd Jx;
    MatrixXd Ju;
    MatrixXd JxN;
}qp_workspace;

typedef struct{
    MatrixXd dx;
    MatrixXd du;
    MatrixXd lam;
    VectorXd mu_u;
    VectorXd mu_x;
    VectorXd mu_g;
}qp_out;

void qp_in_init(const model_size& size, qp_in& in);
void qp_problem_init(const model_size& size, qp_problem& qp);
void qp_workspace_init(const model_size& size, qp_workspace& work);
void qp_generation(const model_size& size, const qp_in& in, qp_workspace& work, qp_problem& qp);
void qp_out_init(model_size& size, qp_out& out);
void expand(model_size& size, qp_in& in, qp_problem& qp, qp_out& out, const VectorXd& x0);

#endif