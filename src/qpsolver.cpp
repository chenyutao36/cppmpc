#include<iostream>
#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include "qpsolver.hpp"
#include <qpOASES.hpp>

USING_NAMESPACE_QPOASES

using namespace std;

qpoases_workspace& qpoases_workspace::init(model_size& size)
{
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    myQP = new SQProblem(N*nu, N*nbx+N*nbg+nbgN);
    myOptions = new Options();
    myOptions->setToMPC();
    myOptions->printLevel = PL_NONE;
    myQP->setOptions(*myOptions);

    return *this;
}

void qpoases_workspace::solveQP(model_size& size, full_condensing_workspace& cond_work,
    qp_problem& qp_data, qp_out& out, int sample)
{
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    int nWSR = 50;

    if(sample==0)
        myQP->init(cond_work.Hc.data(), cond_work.gc.data(), cond_work.Cc.data(), qp_data.lb_u.data(), qp_data.ub_u.data(), cond_work.lcc.data(), cond_work.ucc.data(), nWSR, 0);
    else
        myQP->hotstart(cond_work.Hc.data(), cond_work.gc.data(), cond_work.Cc.data(), qp_data.lb_u.data(), qp_data.ub_u.data(), cond_work.lcc.data(), cond_work.ucc.data(), nWSR, 0);


    double yOpt[N*nu+N*nbx+N*nbg+nbgN];
    myQP->getPrimalSolution(out.du.data());
    myQP->getDualSolution(yOpt);
    memcpy(out.mu_u.data(),yOpt,N*nu*sizeof(double));
    memcpy(out.mu_g.data(),yOpt+N*nu,(N*nbg+nbgN)*sizeof(double));
    memcpy(out.mu_x.data(),yOpt+N*nu+N*nbg+nbgN,N*nbx*sizeof(double));
}

void qpoases_workspace::free()
{
    delete myQP;
    delete myOptions;
}