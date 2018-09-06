#include<iostream>
#include "mpc_common.hpp"
#include "qp_problem.hpp"
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

void qpoases_workspace::solveQP(model_size& size, qp_problem& qp, full_condensing_workspace& cond_work, int sample)
{
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    int nWSR = 50;

    if(sample==0)
        myQP->init(cond_work.Hc.data(), cond_work.gc.data(), cond_work.Cc.data(), qp.data.lb_u.data(), qp.data.ub_u.data(), cond_work.lcc.data(), cond_work.ucc.data(), nWSR, 0);
    else
        myQP->hotstart(cond_work.Hc.data(), cond_work.gc.data(), cond_work.Cc.data(), qp.data.lb_u.data(), qp.data.ub_u.data(), cond_work.lcc.data(), cond_work.ucc.data(), nWSR, 0);


    double yOpt[N*nu+N*nbx+N*nbg+nbgN];
    myQP->getPrimalSolution(qp.out.du.data());
    myQP->getDualSolution(yOpt);
    memcpy(qp.out.mu_u.data(),yOpt,N*nu*sizeof(double));
    memcpy(qp.out.mu_g.data(),yOpt+N*nu,(N*nbg+nbgN)*sizeof(double));
    memcpy(qp.out.mu_x.data(),yOpt+N*nu+N*nbg+nbgN,N*nbx*sizeof(double));

}

qpoases_workspace& qpoases_workspace::free()
{
    delete myQP;
    delete myOptions;

    return *this;
}