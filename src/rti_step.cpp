#include "mpc_common.hpp"
#include "qp_problem.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include "rti_step.hpp"
#include "Timer.h"
#include <iostream>

rti_step_workspace& rti_step_workspace::init(model_size& usr_size)
{
    size = usr_size;
    QP.init(size);
    cond_work.init(size);
    qpoases_work.init(size);
    x0 = VectorXd::Zero(size.nx);
    return *this;
}


rti_step_workspace& rti_step_workspace::step()
{
    timer.start();
    QP.generateQP(size);

    cond_work.full_condensing(size, QP, x0);

    qpoases_work.solveQP(size, QP, cond_work, sample);

    QP.expandSol(size, x0);

    CPT = timer.getElapsedTimeInMilliSec();

    return *this;
}

rti_step_workspace& rti_step_workspace::free()
{
    qpoases_work.free();

    return *this;
}

rti_step_workspace& rti_step_workspace::info()
{
    OBJ = 0;
    QP.info(size,OBJ);
    return *this;
}
