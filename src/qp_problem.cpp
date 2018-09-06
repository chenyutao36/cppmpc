
#include <iostream>
#include <Eigen/Dense>
#include "mpc_common.hpp"
#include "qp_problem.hpp"
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

qp_data& qp_data::init(model_size& size)
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

qp_problem& qp_problem::init(model_size& size)
{
    in.init(size);
    data.init(size);
    work.init(size);
    out.init(size);

    return *this;
}

qp_problem& qp_problem::generateQP(model_size& size)
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
               
    int i=0,j=0;

    // allocate array of pointers   
    double *casadi_in[5];
    double *casadi_out[2];    
    casadi_in[4] = in.W.data();
                      
    // start loop
    for(i=0;i<N;i++){
        casadi_in[0] = in.x.data()+i*nx;
        casadi_in[1] = in.u.data()+i*nu;
        casadi_in[2] = in.p.data()+i*np;
        casadi_in[3] = in.y.data()+i*ny;
        
        // control bounds
        data.lb_u.segment(i*nu,nu) = in.lbu-in.u.col(i);
        data.ub_u.segment(i*nu,nu) = in.ubu-in.u.col(i);
          
        // state bounds
        for (j=0;j<nbx;j++){
            data.lb_x(i*nbx+j) = in.lbx(j)-in.x(nbx_idx[j],i+1);
            data.ub_x(i*nbx+j) = in.ubx(j)-in.x(nbx_idx[j],i+1);
        }
        
        // integration                             
        casadi_out[0] = data.a.data()+i*nx;
        F_Fun(casadi_in, casadi_out);
        // equality residual        
        data.a.col(i) -= in.x.col(i+1);
      
        // sensitivity computation
        casadi_out[0] = data.A.data() + i*nx*nx;
        casadi_out[1] = data.B.data() + i*nx*nu;
        D_Fun(casadi_in, casadi_out);	
                                
        // Hessian
        casadi_out[0] = work.Jx.data();
        casadi_out[1] = work.Ju.data();
        Ji_Fun(casadi_in, casadi_out);

        data.Q.block(0,i*nx,nx,nx) = work.Jx.transpose() * work.Jx;
        data.S.block(0,i*nu,nx,nu) = work.Jx.transpose() * work.Ju;
        data.R.block(0,i*nu,nu,nu) = work.Ju.transpose() * work.Ju;
        regularization(nx, data.Q.data()+i*nx*nx, in.reg);
        regularization(nu, data.R.data()+i*nu*nu, in.reg);
                
        // gradient
        casadi_out[0] = data.gx.data()+i*nx;
        casadi_out[1] = data.gu.data()+i*nu;
        gi_Fun(casadi_in, casadi_out);

        //constraints          
        if (nbg>0){       
            casadi_out[0] = data.lb_g.data() + i*nbg;
            path_con_Fun(casadi_in, casadi_out);
            // constraint residual
            data.ub_g.segment(i*nbg,nbg) = in.ubg - data.lb_g.segment(i*nbg,nbg);
            data.lb_g.segment(i*nbg,nbg) = in.lbg - data.lb_g.segment(i*nbg,nbg);
            
            // constraint Jacobian
            casadi_out[0] = data.Cgx.data()+i*nbg*nx;
            casadi_out[1] = data.Cgu.data()+i*nbg*nu;
            Ci_Fun(casadi_in, casadi_out);
        }
    }
    
    // the terminal stage
    casadi_in[0] = in.x.data()+N*nx;
    casadi_in[1] = in.p.data()+N*np;
    casadi_in[2] = in.yN.data();
    casadi_in[3] = in.WN.data();
    
    casadi_out[0] = work.JxN.data();
    JN_Fun(casadi_in, casadi_out);
    data.Q.block(0,N*nx,nx,nx) = work.JxN.transpose() * work.JxN;
    
    casadi_out[0] = data.gx.data()+N*nx;
    gN_Fun(casadi_in, casadi_out);

    if (nbgN>0){
        casadi_out[0] = data.lb_g.data() + N*nbg;
        path_con_N_Fun(casadi_in, casadi_out);
        data.ub_g.segment(N*nbg,nbgN) = in.ubgN - data.lb_g.segment(N*nbg,nbgN);
        data.lb_g.segment(N*nbg,nbgN) = in.lbgN - data.lb_g.segment(N*nbg,nbgN);

        casadi_out[0] = data.CgN.data();
        CN_Fun(casadi_in, casadi_out);
    }

    return *this;
    
}

qp_problem& qp_problem::expandSol(model_size& size, VectorXd& x0)
 {
    int nx = size.nx;
    int nu = size.nu;
    int N = size.N;

    int i;
    out.dx.col(0) = x0 - in.x.col(0);

    for(i=0;i<N;i++){       
        out.dx.col(i+1) = data.A.block(0,i*nx,nx,nx)*out.dx.col(i)+data.B.block(0,i*nu,nx,nu)*out.du.col(i)+data.a.col(i);
    }

    in.x += out.dx;
    in.u += out.du;

    return *this;

 }

