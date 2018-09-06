#ifndef QPSOLVER_H_
#define QPSOLVER_H_

#include "mpc_common.hpp"
#include "qp_problem.hpp"
#include "full_condensing.hpp"
#include <qpOASES.hpp>

USING_NAMESPACE_QPOASES

class qpoases_workspace{
    public:
        SQProblem *myQP;
        Options *myOptions;

        qpoases_workspace& init(model_size& size);

        void solveQP(model_size& size, qp_problem& qp, full_condensing_workspace& cond_work, int sample);

        qpoases_workspace& free();
};

#endif 