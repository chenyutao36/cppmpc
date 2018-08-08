#include <iostream>
#include <math.h>
#include "mpc_common.hpp"
#include "rti_step.hpp"

using namespace std;
int main()
{  

    // define problem size
    model_size size;
    size.nx = 4;
    size.nu = 1;
    size.ny = 5;
    size.nyN= 4;
    size.np = 0;
    size.nbu = 1;
    size.nbx = 1;
    size.nbg = 0;
    size.nbgN = 0;
    size.N = 5;
    size.nbx_idx=(int*)calloc(size.nbx, sizeof(int));
    size.nbx_idx[0]=0;
    
    // initialize MPC input
    rti_step_workspace rti_work;

    rti_step_init(size, rti_work);

    rti_work.x0(1) = M_PI;

    int i;
    for(i=0;i<size.N+1;i++){
        rti_work.in.x(1,i) = M_PI;
    }
    rti_work.in.W(0,0) = 10; 
    rti_work.in.W(1,1) = 10; 
    rti_work.in.W(2,2) = 0.1; 
    rti_work.in.W(3,3) = 0.1; 
    rti_work.in.W(4,4) = 0.01;

    rti_work.in.WN(0,0) = 10; 
    rti_work.in.WN(1,1) = 10; 
    rti_work.in.WN(2,2) = 0.1; 
    rti_work.in.WN(3,3) = 0.1;

    rti_work.in.lbu(0) = -20; 
    rti_work.in.ubu(0) = 20;
    rti_work.in.lbx(0) = -2; 
    rti_work.in.ubx(0) = 2;

    rti_work.in.reg = 1E-8;

    rti_step(size, rti_work);

    // cout <<"a="<<endl<<rti_work.qp.a<<endl;
    // cout <<"A="<<endl<<rti_work.qp.A.block(0,0,size.nx,size.nx)<<endl;
    // cout <<"B="<<endl<<rti_work.qp.B.block(0,0,size.nx,size.nu)<<endl;
    // cout <<"Q="<<endl<<rti_work.qp.Q.block(0,0,size.nx,size.nx)<<endl;
    // cout <<"S="<<endl<<rti_work.qp.S.block(0,0,size.nx,size.nu)<<endl;
    // cout <<"R="<<endl<<rti_work.qp.R.block(0,0,size.nu,size.nu)<<endl;
    // cout <<"gx="<<endl<<rti_work.qp.gx.col(0)<<endl;
    // cout <<"gu="<<endl<<rti_work.qp.gu.col(0)<<endl;
    // cout <<"lb_u="<<endl<<rti_work.qp.lb_u<<endl;
    // cout <<"lb_x="<<endl<<rti_work.qp.lb_x<<endl;

    // cout <<"Hc="<<endl<<rti_work.cond_work.Hc<<endl;
    // cout <<"Cc="<<endl<<rti_work.cond_work.Cc<<endl;
    // cout <<"gc="<<endl<<rti_work.cond_work.gc<<endl;
    // cout <<"lcc="<<endl<<rti_work.cond_work.lcc<<endl;
    // cout <<"ucc="<<endl<<rti_work.cond_work.ucc<<endl;

    cout <<"xOpt="<<endl<<rti_work.out.dx<<endl;
    cout <<"lam_Opt="<<endl<<rti_work.out.lam<<endl;

    free(size.nbx_idx);

    return 0;
}