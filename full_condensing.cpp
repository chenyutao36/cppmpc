#include "mpc_common.hpp"
#include "qp_generation.hpp"
#include "full_condensing.hpp"

void full_condensing_workspace_init(model_size& size, full_condensing_workspace& cond_work)
{
    int nx = size.nx;
    int nu = size.nu;
    int nbx = size.nbx;
    int nbg = size.nbg;
    int nbgN = size.nbgN;
    int N = size.N;

    cond_work.Hc = Matrix<double, Dynamic, Dynamic, RowMajor>::Zero(N*nu,N*nu);
    cond_work.Cc = Matrix<double, Dynamic, Dynamic, RowMajor>::Zero(N*nbx+N*nbg+nbgN,N*nu);
    cond_work.gc = VectorXd::Zero(N*nu);
    cond_work.lcc = VectorXd::Zero(N*nbx+N*nbg+nbgN);
    cond_work.ucc = VectorXd::Zero(N*nbx+N*nbg+nbgN);

    cond_work.G = MatrixXd::Zero(N*nx,N*nu);
    cond_work.L = VectorXd::Zero((N+1)*nx);
    cond_work.W = MatrixXd::Zero(nx,nu);
    cond_work.w = VectorXd::Zero(nx);
}

void full_condensing(model_size& size, full_condensing_workspace& cond_work,
    qp_in& in, qp_problem& qp, const VectorXd& x0)
{
    int nx = size.nx;
    int nu = size.nu;
    int nbx = size.nbx;
    int nbg = size.nbg;
    int nbgN = size.nbgN;
    int N = size.N;

    int idx = N*nbg+nbgN;      

    int i,j;

    /* compute G */
    for(i=0;i<N;i++){
        cond_work.G.block(i*nx,i*nu,nx,nu) = qp.B.block(0,i*nu,nx,nu);
        for(j=i+1;j<N;j++)
            cond_work.G.block(j*nx,i*nu,nx,nu) = qp.A.block(0,j*nx,nx,nx) * cond_work.G.block((j-1)*nx,i*nu,nx,nu);
    }

    /* Compute Hc */
    for(i=0;i<N;i++){
        cond_work.W = qp.Q.block(0,N*nx,nx,nx) * cond_work.G.block((N-1)*nx,i*nu,nx,nu);
        for(j=N-1;j>i;j--){     
            cond_work.Hc.block(j*nu,i*nu,nu,nu) = qp.S.block(0,j*nu,nx,nu).transpose() * cond_work.G.block((j-1)*nx,i*nu,nx,nu) + qp.B.block(0,j*nu,nx,nu).transpose() * cond_work.W; 
            cond_work.W = qp.Q.block(0,j*nx,nx,nx) * cond_work.G.block((j-1)*nx,i*nu,nx,nu) + qp.A.block(0,j*nx,nx,nx).transpose() * cond_work.W;
            cond_work.Hc.block(i*nu,j*nu,nu,nu) = cond_work.Hc.block(j*nu,i*nu,nu,nu).transpose();
        }
        cond_work.Hc.block(i*nu,i*nu,nu,nu) = qp.R.block(0,i*nu,nu,nu) + qp.B.block(0,i*nu,nx,nu).transpose() * cond_work.W;
    }

    /* Compute Ccg */
    if (nbg>0){         
        for(i=0;i<N;i++){
            cond_work.Cc.block(i*nbg,i*nu,nbg,nu) = qp.Cgu.block(0,i*nu,nbg,nu);
            for(j=i+1;j<N;j++){   
                cond_work.Cc.block(j*nbg,i*nu,nbg,nu) = qp.Cgx.block(0,j*nx,nbg,nx)*cond_work.G.block((j-1)*nx,i*nu,nx,nu);
            }    
        }  
    }
    /* Compute CcN */
    if (nbgN>0){          
        for(i=0;i<N;i++){                 
            cond_work.Cc.block(N*nbg,i*nu,nbgN,nu) = qp.CgN * cond_work.G.block((N-1)*nx,i*nu,nx,nu);
        }
    }

    /* Compute Ccx */
    if (nbx>0){   
        for(i=0;i<N;i++){
            for(j=i+1;j<=N;j++) 
                cond_work.Cc.block(idx+(j-1)*nbx,i*nu,nbx,nu) = qp.Cx * cond_work.G.block((j-1)*nx,i*nu,nx,nu);
                
        }  
    }

    /* compute L */
    cond_work.L.head(nx) = x0 - in.x.col(0);
    for(i=0;i<N;i++){
        cond_work.L.segment((i+1)*nx,nx) = qp.A.block(0,i*nx,nx,nx)*cond_work.L.segment(i*nx,nx) + qp.a.col(i);
    }

    /* compute gc */
    cond_work.w = qp.gx.col(N) + qp.Q.block(0,N*nx,nx,nx) * cond_work.L.tail(nx);
    for(i=N-1;i>0;i--){
        cond_work.gc.segment(i*nu,nu) = qp.gu.col(i) + qp.S.block(0,i*nu,nx,nu).transpose()*cond_work.L.segment(i*nx,nx) + qp.B.block(0,i*nu,nx,nu).transpose()*cond_work.w;
        cond_work.w = qp.gx.col(i) + qp.Q.block(0,i*nx,nx,nx)*cond_work.L.segment(i*nx,nx) + qp.A.block(0,i*nx,nx,nx).transpose()*cond_work.w;
    }
    cond_work.gc.head(nu) = qp.gu.col(0) + qp.S.leftCols(nu).transpose()*cond_work.L.head(nx) + qp.B.leftCols(nu).transpose()*cond_work.w;

    /* Compute cc */
    if (nbg>0){ 
        for(i=0;i<N;i++){
            cond_work.lcc.segment(i*nbg,nbg) = qp.Cgx.block(0,i*nx,nbg,nx)*cond_work.L.segment(i*nx,nx);
            cond_work.ucc.segment(i*nbg,nbg) = qp.ub_g.segment(i*nbg,nbg) - cond_work.lcc.segment(i*nbg,nbg);
            cond_work.lcc.segment(i*nbg,nbg) = qp.lb_g.segment(i*nbg,nbg) - cond_work.lcc.segment(i*nbg,nbg);
        }       
    }
    /* Compute ccN */
    if (nbgN>0){ 
        cond_work.lcc.segment(N*nbg,nbgN) = qp.CgN * cond_work.L.tail(nx);
        cond_work.ucc.segment(N*nbg,nbgN) = qp.ub_g.segment(N*nbg,nbgN) - cond_work.lcc.segment(N*nbg,nbgN);
        cond_work.lcc.segment(N*nbg,nbgN) = qp.lb_g.segment(N*nbg,nbgN) - cond_work.lcc.segment(N*nbg,nbgN);
    }

    if (nbx>0){       
        for(i=0;i<N;i++){
            cond_work.lcc.segment(idx+i*nbx,nbx) = qp.Cx * cond_work.L.segment((i+1)*nx,nx);
            cond_work.ucc.segment(idx+i*nbx,nbx) = qp.ub_x.segment(i*nbx,nbx) - cond_work.lcc.segment(idx+i*nbx,nbx);
            cond_work.lcc.segment(idx+i*nbx,nbx) = qp.lb_x.segment(i*nbx,nbx) - cond_work.lcc.segment(idx+i*nbx,nbx);
        }
    }
}