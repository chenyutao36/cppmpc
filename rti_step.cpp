#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include "rti_step.hpp"

void rti_step_init(model_size& size, rti_step_workspace& rti_work)
{
    qp_in_init(size, rti_work.in);
    qp_problem_init(size, rti_work.qp);
    qp_workspace_init(size, rti_work.qp_work);
    full_condensing_workspace_init(size, rti_work.cond_work);
    qp_out_init(size, rti_work.out);
    qpoases_workspace_init(size,rti_work.qpoases_work);
    rti_work.x0 = VectorXd::Zero(size.nx);
}

void rti_step(model_size& size, rti_step_workspace& rti_work)
{
    qp_generation(size, rti_work.in, rti_work.qp_work, rti_work.qp);

    full_condensing(size, rti_work.cond_work, rti_work.in, rti_work.qp, rti_work.x0);

    solveQP(size, rti_work.cond_work, rti_work.qp, rti_work.qpoases_work, rti_work.out,
     rti_work.sample);

    expand(size, rti_work.in, rti_work.qp, rti_work.out, rti_work.x0);

}

void rti_step_free(rti_step_workspace& rti_work)
{
    qpoases_workspace_free(rti_work.qpoases_work);
}