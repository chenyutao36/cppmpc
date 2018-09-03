
#include <iostream>
#include <Eigen/Dense>
#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "casadi_wrapper.hpp"

using namespace Eigen;
using namespace std;

qp_in& qp_in::init(model_size& size)
{
    int nx = size.nx;
    int nu = size.nu;
    int ny = size.ny;
    int nyN = size.nyN;
    int np = size.np;
    int nbx = size.nbx;
    int nbu = size.nbu;
    int nbg = size.nbg;
    int nbgN = size.nbgN;
    int N = size.N;
    
    x = MatrixXd::Zero(nx,N+1);
    u = MatrixXd::Zero(nu,N);
    y = MatrixXd::Zero(ny,N);
    yN = VectorXd::Zero(nyN);
    p = MatrixXd::Zero(np,N+1);
    W = MatrixXd::Zero(ny,ny);
    WN = MatrixXd::Zero(nyN,nyN);
    lbu = VectorXd::Zero(nbu);
    ubu = VectorXd::Zero(nbu);
    lbx = VectorXd::Zero(nbx);
    ubx = VectorXd::Zero(nbx);
    lbg = VectorXd::Zero(nbg);
    ubg = VectorXd::Zero(nbg);
    lbgN = VectorXd::Zero(nbgN);
    ubgN = VectorXd::Zero(nbgN);

    return *this;
}

qp_problem& qp_problem::init(model_size& size)
{
    int nx=size.nx;
    int nu=size.nu;
    int nbx = size.nbx;
    int nbg = size.nbg;
    int nbgN = size.nbgN;
    int N = size.N;
    int *nbx_idx = size.nbx_idx;

    Q = MatrixXd::Zero(nx,(N+1)*nx);
    S = MatrixXd::Zero(nx,N*nu);
    R = MatrixXd::Zero(nu,N*nu);
    A = MatrixXd::Zero(nx,N*nx);
    B = MatrixXd::Zero(nx,N*nu);
    a = MatrixXd::Zero(nx,N);
    Cx = MatrixXd::Zero(nbx,nx);
    Cgx = MatrixXd::Zero(nbg,N*nx);
    CgN = MatrixXd::Zero(nbgN,nx);
    Cgu = MatrixXd::Zero(nbg,N*nu);
    gx = MatrixXd::Zero(nx,N+1);
    gu = MatrixXd::Zero(nu,N);
    lb_u = VectorXd::Zero(N*nu);
    ub_u = VectorXd::Zero(N*nu);
    lb_x = VectorXd::Zero(N*nbx);
    ub_x = VectorXd::Zero(N*nbx);
    lb_g = VectorXd::Zero(nbg*N+nbgN);
    ub_g = VectorXd::Zero(nbg*N+nbgN);

    int i;
    for(i=0;i<nbx;i++)
        Cx(i,nbx_idx[i]) = 1.0;

    return *this;

}

qp_workspace& qp_workspace::init(model_size& size)
{
    int nx=size.nx;
    int nu=size.nu;
    int ny=size.ny;
    int nyN=size.nyN;

    Jx = MatrixXd::Zero(ny,nx);
    Ju = MatrixXd::Zero(ny,nu);
    JxN = MatrixXd::Zero(nyN,nx);

    return *this;

}

qp_out& qp_out::init(model_size& size)
{
    int nx=size.nx;
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    dx = MatrixXd::Zero(nx,N+1);
    du = MatrixXd::Zero(nu,N);
    lam = MatrixXd::Zero(nx,N+1);
    mu_u = VectorXd::Zero(N*nu);   
    mu_x = VectorXd::Zero(N*nbx);
    mu_g = VectorXd::Zero(nbg*N+nbgN);

    return *this;

}

void qp_generation(model_size& size, qp_in& in, qp_workspace& work, qp_problem& qp)
{
    int nx = size.nx;
    int nu = size.nu;
    int ny = size.ny;
    int np = size.np;
    int nbx = size.nbx;
    int nbg = size.nbg;
    int nbgN = size.nbgN;
    int N = size.N;
    int *nbx_idx = size.nbx_idx;

    MatrixXd x = in.x;
    MatrixXd u = in.u;
    MatrixXd y = in.y;
    VectorXd yN = in.yN; 
    MatrixXd W = in.W;
    MatrixXd WN = in.WN;
    MatrixXd p = in.p;
    VectorXd lbu = in.lbu;
    VectorXd ubu = in.ubu;
    VectorXd lbx = in.lbx;
    VectorXd ubx = in.ubx;
    VectorXd lbg = in.lbg;
    VectorXd ubg = in.ubg;
    VectorXd lbgN = in.lbgN;
    VectorXd ubgN = in.ubgN;
    double reg = in.reg;
               
    int i=0,j=0;

    // allocate array of pointers   
    double *casadi_in[5];
    double *casadi_out[2];    
    casadi_in[4] = W.data();
                      
    // start loop
    for(i=0;i<N;i++){
        casadi_in[0] = x.data()+i*nx;
        casadi_in[1] = u.data()+i*nu;
        casadi_in[2] = p.data()+i*np;
        casadi_in[3] = y.data()+i*ny;
        
        // control bounds
        qp.lb_u.segment(i*nu,nu) = lbu-u.col(i);
        qp.ub_u.segment(i*nu,nu) = ubu-u.col(i);
          
        // state bounds
        for (j=0;j<nbx;j++){
            qp.lb_x(i*nbx+j) = lbx(j)-x(nbx_idx[j],i+1);
            qp.ub_x(i*nbx+j) = ubx(j)-x(nbx_idx[j],i+1);
        }
        
        // integration                             
        casadi_out[0] = qp.a.data()+i*nx;
        F_Fun(casadi_in, casadi_out);
        // equality residual        
        qp.a.col(i) -= x.col(i+1);
      
        // sensitivity computation
        casadi_out[0] = qp.A.data() + i*nx*nx;
        casadi_out[1] = qp.B.data() + i*nx*nu;
        D_Fun(casadi_in, casadi_out);	
                                
        // Hessian
        casadi_out[0] = work.Jx.data();
        casadi_out[1] = work.Ju.data();
        Ji_Fun(casadi_in, casadi_out);

        qp.Q.block(0,i*nx,nx,nx) = work.Jx.transpose() * work.Jx;
        qp.S.block(0,i*nu,nx,nu) = work.Jx.transpose() * work.Ju;
        qp.R.block(0,i*nu,nu,nu) = work.Ju.transpose() * work.Ju;
        regularization(nx, qp.Q.data()+i*nx*nx, reg);
        regularization(nu, qp.R.data()+i*nu*nu, reg);
                
        // gradient
        casadi_out[0] = qp.gx.data()+i*nx;
        casadi_out[1] = qp.gu.data()+i*nu;
        gi_Fun(casadi_in, casadi_out);

        //constraints          
        if (nbg>0){       
            casadi_out[0] = qp.lb_g.data() + i*nbg;
            path_con_Fun(casadi_in, casadi_out);
            // constraint residual
            qp.ub_g.segment(i*nbg,nbg) = ubg - qp.lb_g.segment(i*nbg,nbg);
            qp.lb_g.segment(i*nbg,nbg) = lbg - qp.lb_g.segment(i*nbg,nbg);
            
            // constraint Jacobian
            casadi_out[0] = qp.Cgx.data()+i*nbg*nx;
            casadi_out[1] = qp.Cgu.data()+i*nbg*nu;
            Ci_Fun(casadi_in, casadi_out);
        }
    }
    
    // the terminal stage
    casadi_in[0] = x.data()+N*nx;
    casadi_in[1] = p.data()+N*np;
    casadi_in[2] = yN.data();
    casadi_in[3] = WN.data();
    
    casadi_out[0] = work.JxN.data();
    JN_Fun(casadi_in, casadi_out);
    qp.Q.block(0,N*nx,nx,nx) = work.JxN.transpose() * work.JxN;
    
    casadi_out[0] = qp.gx.data()+N*nx;
    gN_Fun(casadi_in, casadi_out);

    if (nbgN>0){
        casadi_out[0] = qp.lb_g.data() + N*nbg;
        path_con_N_Fun(casadi_in, casadi_out);
        qp.ub_g.segment(N*nbg,nbgN) = ubgN - qp.lb_g.segment(N*nbg,nbgN);
        qp.lb_g.segment(N*nbg,nbgN) = lbgN - qp.lb_g.segment(N*nbg,nbgN);

        casadi_out[0] = qp.CgN.data();
        CN_Fun(casadi_in, casadi_out);
    }
    
}

void expand(model_size& size, qp_in& in, qp_problem& qp, qp_out& out, VectorXd& x0)
 {
    int nx = size.nx;
    int nu = size.nu;
    int N = size.N;

    int i;
    out.dx.col(0) = x0 - in.x.col(0);

    for(i=0;i<N;i++){       
        out.dx.col(i+1) = qp.A.block(0,i*nx,nx,nx)*out.dx.col(i)+qp.B.block(0,i*nu,nx,nu)*out.du.col(i)+qp.a.col(i);
    }

    in.x += out.dx;
    in.u += out.du;

 }

