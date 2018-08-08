#ifndef RTI_STEP_H_
#define RTI_STEP_H_

#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"

typedef struct{
    qp_in in;
    qp_problem qp;
    qp_workspace qp_work;  
    qp_out out;

    full_condensing_workspace cond_work;

    VectorXd x0;
}rti_step_workspace;

void rti_step_init(model_size& size, rti_step_workspace& rti_work);

void rti_step(model_size& size, rti_step_workspace& rti_work);

#endif