#include <iostream>
#include <fstream>
#include <math.h>
#include "mpc_common.hpp"
#include "rti_step.hpp"
#include "casadi_wrapper.hpp"

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
    size.N = 40;
    size.nbx_idx = new int[size.nbx];            
    size.nbx_idx[0]=0;
    
    // create a RTI controller
    rti_step_workspace rti_work(size);

    // initial condition and parameters
    rti_work.x0(1) = M_PI;

    int i;
    for(i=0;i<size.N+1;i++){
        rti_work.QP.in.x(1,i) = M_PI;
    }
    rti_work.QP.in.W(0,0) = 10; 
    rti_work.QP.in.W(1,1) = 10; 
    rti_work.QP.in.W(2,2) = 0.1; 
    rti_work.QP.in.W(3,3) = 0.1; 
    rti_work.QP.in.W(4,4) = 0.01;

    rti_work.QP.in.WN(0,0) = 10; 
    rti_work.QP.in.WN(1,1) = 10; 
    rti_work.QP.in.WN(2,2) = 0.1; 
    rti_work.QP.in.WN(3,3) = 0.1;

    rti_work.QP.in.lbu(0) = -20; 
    rti_work.QP.in.ubu(0) = 20;
    rti_work.QP.in.lbx(0) = -2; 
    rti_work.QP.in.ubx(0) = 2;

    rti_work.QP.in.reg = 1E-8;

    // prepare the closed-loop simulation
    rti_work.sample = 0;
    double Tf=4, Ts=0.05,t=0;

    double *simu_in[3];
    double *simu_out[1];
    simu_in[2] = rti_work.QP.in.p.col(0).data();

    ofstream myfile;
    myfile.open ("data.txt");

    // start the simulation
    while(t<Tf){
        
        // call RTI solving routine
        rti_work.step();

        // feedback
        simu_in[0] = rti_work.QP.in.x.col(0).data();
        simu_in[1] = rti_work.QP.in.u.col(0).data();
        simu_out[0] = rti_work.x0.data();
        F_Fun(simu_in, simu_out);

        // update sampling and time
        rti_work.sample++;
        t += Ts;

        // store the closed-loop results
        myfile <<"Sample " << rti_work.sample <<": " << rti_work.x0.transpose() << " | " << rti_work.QP.in.u.col(0) <<" |CPT=: " << rti_work.CPT << "ms" <<endl;

        // shifting(optional)
        for(i=0;i<size.N-1;i++){
            rti_work.QP.in.x.col(i) = rti_work.QP.in.x.col(i+1);
            rti_work.QP.in.u.col(i) = rti_work.QP.in.u.col(i+1);
        }
        rti_work.QP.in.x.col(size.N-1) = rti_work.QP.in.x.col(size.N);   

    }

    // free memory
    myfile.close();
      
    delete [] size.nbx_idx;
    size.nbx_idx = NULL;
    
    rti_work.free();

    return 0;
}