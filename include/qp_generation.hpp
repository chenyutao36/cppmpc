#ifndef QP_GENERATION_H_
#define QP_GENERATION_H_

#include <Eigen/Dense>
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

class qp_problem{
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
        VectorXd lb_u;
        VectorXd ub_u;
        VectorXd lb_x;
        VectorXd ub_x;
        VectorXd lb_g;
        VectorXd ub_g;
    
        qp_problem& init(model_size& size);
};

class qp_workspace{
    public:
        MatrixXd Jx;
        MatrixXd Ju;
        MatrixXd JxN;

        qp_workspace& init(model_size& size);
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

void qp_generation(model_size& size, qp_in& in, qp_workspace& work, qp_problem& qp);
void expand(model_size& size, qp_in& in, qp_problem& qp, qp_out& out, VectorXd& x0);

#endif