#include<iostream>
#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include <qpOASES.hpp>

USING_NAMESPACE_QPOASES

using namespace std;

void qpoases_workspace_init(const model_size& size, qpoases_workspace &qpoases_work)
{
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    qpoases_work.myQP = new SQProblem(N*nu, N*nbx+N*nbg+nbgN);
    qpoases_work.myOptions = new Options();
    qpoases_work.myOptions->setToMPC();
    qpoases_work.myOptions->printLevel = PL_NONE;
    qpoases_work.myQP->setOptions(*qpoases_work.myOptions);
}

void solveQP(const model_size& size, full_condensing_workspace& cond_work,
    qp_problem& qp_data, qpoases_workspace &qpoases_work, qp_out& out,
    int sample)
{
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    int nWSR = 50;

    if(sample==0)
        qpoases_work.myQP->init(cond_work.Hc.data(), cond_work.gc.data(), cond_work.Cc.data(), qp_data.lb_u.data(), qp_data.ub_u.data(), cond_work.lcc.data(), cond_work.ucc.data(), nWSR, 0);
    else
        qpoases_work.myQP->hotstart(cond_work.Hc.data(), cond_work.gc.data(), cond_work.Cc.data(), qp_data.lb_u.data(), qp_data.ub_u.data(), cond_work.lcc.data(), cond_work.ucc.data(), nWSR, 0);


    double yOpt[N*nu+N*nbx+N*nbg+nbgN];
    qpoases_work.myQP->getPrimalSolution(out.du.data());
    qpoases_work.myQP->getDualSolution(yOpt);
    memcpy(out.mu_u.data(),yOpt,N*nu*sizeof(double));
    memcpy(out.mu_g.data(),yOpt+N*nu,(N*nbg+nbgN)*sizeof(double));
    memcpy(out.mu_x.data(),yOpt+N*nu+N*nbg+nbgN,N*nbx*sizeof(double));
}

void qpoases_workspace_free(qpoases_workspace &qpoases_work)
{
    delete qpoases_work.myQP;
    delete qpoases_work.myOptions;
}