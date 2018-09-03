#ifndef RTI_STEP_H_
#define RTI_STEP_H_

#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include "Timer.h"

class rti_step_workspace{
    public:
        rti_step_workspace() = default;
        rti_step_workspace(model_size& size) {init(size);};

        qp_in in;      
        VectorXd x0;
        int sample;      
        double CPT;
        
        rti_step_workspace& init(model_size& size);
        rti_step_workspace& step(model_size& size);
        rti_step_workspace& free();

    private:
        Timer timer;
        qp_problem qp;
        qp_workspace qp_work; 
        qp_out out;
        full_condensing_workspace cond_work;
        qpoases_workspace qpoases_work; 

};

#endif