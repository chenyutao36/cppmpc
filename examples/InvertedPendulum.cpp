#include <iostream>
#include <fstream>
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
    size.nyN = 4;
    size.np = 0;
    size.nbu = 1;
    size.nbx = 1;
    size.nbg = 0;
    size.nbgN = 0;
    size.N = 40;
    size.nbx_idx = new int[size.nbx];
    size.nbx_idx[0] = 0;
    size.nbu_idx = new int[size.nbu];
    size.nbu_idx[0] = 0;

    // create a RTI controller
    rti_step_workspace rti_work(size);

    // initial condition and parameters
    rti_work.x0(1) = M_PI;

    for(int i=0;i<size.N+1;i++)
        rti_work.QP.in.x(1,i) = M_PI;

    rti_work.QP.in.W(0,0) = 10;
    rti_work.QP.in.W(1,0) = 10;
    rti_work.QP.in.W(2,0) = 0.1;
    rti_work.QP.in.W(3,0) = 0.1;
    rti_work.QP.in.W(4,0) = 0.01;
    for(int i=1;i<size.N;i++)
        rti_work.QP.in.W.col(i) = rti_work.QP.in.W.col(0);

    rti_work.QP.in.WN(0) = 10;
    rti_work.QP.in.WN(1) = 10;
    rti_work.QP.in.WN(2) = 0.1;
    rti_work.QP.in.WN(3) = 0.1;

    rti_work.QP.in.lbu(0) = -20;
    rti_work.QP.in.ubu(0) = 20;
    rti_work.QP.in.lbx(0) = -2;
    rti_work.QP.in.ubx(0) = 2;

    rti_work.QP.in.reg = 1E-8;

    // prepare the closed-loop simulation
    rti_work.sample = 0;
    double Tf=4, Ts=0.025,t=0;

    double *simu_in[3];
    double *simu_out[1];
    simu_in[2] = rti_work.QP.in.p.col(0).data();

    ofstream myfile;
    myfile.open ("data.txt");

    rti_work.info();
    myfile <<"Sample 0: " << rti_work.x0.transpose() << " | - |OBJ=: "<<rti_work.OBJ <<" |CPT=: " << "" <<endl;
    // start the simulation
    while (t < Tf)
    {
        // call RTI solving routine
        rti_work.step();
        rti_work.info();

        // feedback
        simu_in[0] = rti_work.QP.in.x.col(0).data();
        simu_in[1] = rti_work.QP.in.u.col(0).data();
        simu_out[0] = rti_work.x0.data();
        F_Fun(simu_in, simu_out);

        // update sampling and time
        rti_work.sample++;
        t += Ts;

        // store the closed-loop results
        myfile <<"Sample " << rti_work.sample <<": " << rti_work.x0.transpose() << " | " << rti_work.QP.in.u.col(0) << " |OBJ=: "<<rti_work.OBJ <<" |CPT=: " << rti_work.CPT << "ms" <<endl;

        // shifting(optional)
        for (int i=0; i<size.N-1; i++)
        {
            rti_work.QP.in.x.col(i) = rti_work.QP.in.x.col(i+1);
            rti_work.QP.in.u.col(i) = rti_work.QP.in.u.col(i+1);
        }
        rti_work.QP.in.x.col(size.N-1) = rti_work.QP.in.x.col(size.N);
    }

    // free memory
    myfile.close();

    delete [] size.nbx_idx;
    size.nbx_idx = NULL;
    size.nbu_idx = NULL;

    rti_work.free();

    return 0;
}
