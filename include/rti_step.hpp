#ifndef RTI_STEP_H_
#define RTI_STEP_H_

#include "mpc_common.hpp"
#include "qp_problem.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include "Timer.h"

class rti_step_workspace{
    public:
        rti_step_workspace() = default;
        rti_step_workspace(model_size& size) {init(size);};

        model_size size; // problem dimension
        qp_problem QP;   // QP subproblem
        VectorXd x0;     // current state
        int sample;      // current sample
        double CPT;      // computation time
        double OBJ;
        rti_step_workspace& init(model_size& size);
        rti_step_workspace& step();
        rti_step_workspace& free();
        rti_step_workspace& info();

    private:
        Timer timer;
        full_condensing_workspace cond_work;
        qpoases_workspace qpoases_work;
};

#endif
