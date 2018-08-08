#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"
#include "qpOASES.hpp"
#include "qpsolver.hpp"

USING_NAMESPACE_QPOASES

void solveQP(const model_size& size, full_condensing_workspace& cond_work,
    qp_problem& qp_data, qp_out& out)
{
    int nx=size.nx;
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    // qpsolver.myQP
    SQProblem myQP(N*nu, N*nbx+N*nbg+nbgN);
    // Options options;
    // myQP.setOptions( options );

    int nWSR = 20;
    myQP.init(cond_work.Hc.data(), cond_work.gc.data(), cond_work.Cc.data(), qp_data.lb_u.data(), qp_data.ub_u.data(), cond_work.lcc.data(), cond_work.ucc.data(), nWSR, 0);

    double yOpt[N*nu+N*nbx+N*nbg+nbgN];
    myQP.getPrimalSolution(out.dx.data());
    myQP.getPrimalSolution(yOpt);

    std::memcpy(out.lam.data(),yOpt,(N+1)*nx);
}