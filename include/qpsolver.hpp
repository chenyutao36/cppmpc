#ifndef QPSOLVER_H_
#define QPSOLVER_H_

#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include <qpOASES.hpp>

USING_NAMESPACE_QPOASES

typedef struct{
    SQProblem *myQP;
    Options *myOptions;
}qpoases_workspace;

void qpoases_workspace_init(const model_size& size, qpoases_workspace &qpoases_work);

void solveQP(const model_size& size, full_condensing_workspace& cond_work,
    qp_problem& qp_data, qpoases_workspace &qpoases_work, qp_out& out,
    int sample);

void qpoases_workspace_free(qpoases_workspace &qpoases_work);

#endif 