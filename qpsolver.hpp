#ifndef QPSOLVER_H_
#define QPSOLVER_H_

#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"

// USING_NAMESPACE_QPOASES

// typedef struct{
//     QProblem myQP;
// }qpoases_workspace;

void solveQP(const model_size& size, full_condensing_workspace& cond_work,
    qp_problem& qp_data, qp_out& out);

#endif 