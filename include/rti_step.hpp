#ifndef RTI_STEP_H_
#define RTI_STEP_H_

#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include "Timer.h"

typedef struct{
    qp_in in;
    qp_problem qp;
    qp_workspace qp_work; 
    qp_out out;

    full_condensing_workspace cond_work;
    qpoases_workspace qpoases_work; 

    VectorXd x0;

    int sample;
    Timer timer;
    double CPT;
}rti_step_workspace;

void rti_step_init(model_size& size, rti_step_workspace& rti_work);

void rti_step(model_size& size, rti_step_workspace& rti_work);

void rti_step_free(rti_step_workspace& rti_work);

#endif