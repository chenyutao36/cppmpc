#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include "rti_step.hpp"
#include "Timer.h"

rti_step_workspace& rti_step_workspace::init(model_size& size)
{
    in.init(size);
    qp.init(size);
    qp_work.init(size);
    cond_work.init(size);
    out.init(size);
    qpoases_work.init(size);
    x0 = VectorXd::Zero(size.nx);

    return *this;
}

rti_step_workspace& rti_step_workspace::step(model_size& size)
{
    timer.start();    

    qp_generation(size, in, qp_work, qp);

    cond_work.full_condensing(size, in, qp, x0);

    qpoases_work.solveQP(size, cond_work, qp, out, sample);

    expand(size, in, qp, out, x0);

    timer.stop();
    CPT = timer.getElapsedTimeInMilliSec();  

    return *this; 
}

rti_step_workspace& rti_step_workspace::free()
{
    qpoases_work.free();

    return *this;
}