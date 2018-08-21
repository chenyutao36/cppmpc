
#include <iostream>
#include <Eigen/Dense>
#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "casadi_wrapper.hpp"

using namespace Eigen;
using namespace std;

void qp_in_init(const model_size& size, qp_in& in)
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
    
    in.x = MatrixXd::Zero(nx,N+1);
    in.u = MatrixXd::Zero(nu,N);
    in.y = MatrixXd::Zero(ny,N);
    in.yN = VectorXd::Zero(nyN);
    in.p = MatrixXd::Zero(np,N+1);
    in.W = MatrixXd::Zero(ny,ny);
    in.WN = MatrixXd::Zero(nyN,nyN);
    in.lbu = VectorXd::Zero(nbu);
    in.ubu = VectorXd::Zero(nbu);
    in.lbx = VectorXd::Zero(nbx);
    in.ubx = VectorXd::Zero(nbx);
    in.lbg = VectorXd::Zero(nbg);
    in.ubg = VectorXd::Zero(nbg);
    in.lbgN = VectorXd::Zero(nbgN);
    in.ubgN = VectorXd::Zero(nbgN);
}



void qp_problem_init(const model_size& size, qp_problem& qp)
{
    int nx=size.nx;
    int nu=size.nu;
    int nbx = size.nbx;
    int nbu = size.nbu;
    int nbg = size.nbg;
    int nbgN = size.nbgN;
    int N = size.N;
    int *nbx_idx = size.nbx_idx;

    qp.Q = MatrixXd::Zero(nx,(N+1)*nx);
    qp.S = MatrixXd::Zero(nx,N*nu);
    qp.R = MatrixXd::Zero(nu,N*nu);
    qp.A = MatrixXd::Zero(nx,N*nx);
    qp.B = MatrixXd::Zero(nx,N*nu);
    qp.a = MatrixXd::Zero(nx,N);
    qp.Cx = MatrixXd::Zero(nbx,nx);
    qp.Cgx = MatrixXd::Zero(nbg,N*nx);
    qp.CgN = MatrixXd::Zero(nbgN,nx);
    qp.Cgu = MatrixXd::Zero(nbg,N*nu);
    qp.gx = MatrixXd::Zero(nx,N+1);
    qp.gu = MatrixXd::Zero(nu,N);
    qp.lb_u = VectorXd::Zero(N*nu);
    qp.ub_u = VectorXd::Zero(N*nu);
    qp.lb_x = VectorXd::Zero(N*nbx);
    qp.ub_x = VectorXd::Zero(N*nbx);
    qp.lb_g = VectorXd::Zero(nbg*N+nbgN);
    qp.ub_g = VectorXd::Zero(nbg*N+nbgN);

    int i;
    for(i=0;i<nbx;i++)
        qp.Cx(i,nbx_idx[i]) = 1.0;

}

void qp_workspace_init(const model_size& size, qp_workspace& work)
{
    int nx=size.nx;
    int nu=size.nu;
    int ny=size.ny;
    int nyN=size.nyN;

    work.Jx = MatrixXd::Zero(ny,nx);
    work.Ju = MatrixXd::Zero(ny,nu);
    work.JxN = MatrixXd::Zero(nyN,nx);

}

void qp_out_init(model_size& size, qp_out& out)
{
    int nx=size.nx;
    int nu=size.nu; 
    int nbg=size.nbg;
    int nbgN=size.nbgN;  
    int N = size.N;
    int nbx = size.nbx;

    out.dx = MatrixXd::Zero(nx,N+1);
    out.du = MatrixXd::Zero(nu,N);
    out.lam = MatrixXd::Zero(nx,N+1);
    out.mu_u = VectorXd::Zero(N*nu);   
    out.mu_x = VectorXd::Zero(N*nbx);
    out.mu_g = VectorXd::Zero(nbg*N+nbgN);

}

void qp_generation(const model_size& size, const qp_in& in, qp_workspace& work, qp_problem& qp)
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
    int *nbx_idx = size.nbx_idx;
    int *nbu_idx = size.nbu_idx;

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

void expand(model_size& size, qp_in& in, qp_problem& qp, qp_out& out, const VectorXd& x0)
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

