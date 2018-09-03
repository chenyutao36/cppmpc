#ifndef FULL_CONDENSING_H_
#define FULL_CONDENSING_H_

#include "mpc_common.hpp"
#include "qp_generation.hpp"

using namespace Eigen;

class full_condensing_workspace{
    public: 
        Matrix<double, Dynamic, Dynamic, RowMajor> Hc;
        Matrix<double, Dynamic, Dynamic, RowMajor> Cc;
        VectorXd gc;
        VectorXd lcc;
        VectorXd ucc;

        MatrixXd G;
        VectorXd L;
        MatrixXd W;
        VectorXd w;

        full_condensing_workspace& init(model_size& size);

        void full_condensing(model_size& size, qp_in& in, qp_problem& qp, VectorXd& x0);
};


#endif