#ifndef QP_PROBLEM_H_
#define QP_PROBLEM_H_

#include <eigen3/Eigen/Dense>
#include "mpc_common.hpp"

using namespace Eigen;

class qp_in{
    public:
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

        qp_in& init(model_size& size);
};

class qp_data{
    public:
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
        VectorXd lb_x;
        VectorXd ub_x;
        VectorXd lb_u;
        VectorXd ub_u;
        VectorXd lb_g;
        VectorXd ub_g;

        qp_data& init(model_size& size);
};

class qp_out{
    public:
        MatrixXd dx;
        MatrixXd du;
        MatrixXd lam;
        VectorXd mu_u;
        VectorXd mu_x;
        VectorXd mu_g;

        qp_out& init(model_size& size);
};

class qp_problem{
    public:
        qp_in in;
        qp_data data;
        qp_out out;

        qp_problem& init(model_size& size);
        qp_problem& generateQP(model_size& size);
        qp_problem& expandSol(model_size& size, VectorXd& x0);
        qp_problem& info(model_size& size, double& OBJ);
};

#endif
