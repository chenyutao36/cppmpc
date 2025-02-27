#include <stdlib.h>
#include "casadi_wrapper.hpp"
#include "casadi_src.h"

void f_Fun(double **in, double **out){
    const double *x = in[0];
    const double *u = in[1];
    const double *p = in[2];

    double *xn = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;

    casadi_res[0] = xn;

    f_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void vde_Fun(double **in, double **out){
    const double *x = in[0];
    const double *u = in[1];
    const double *p = in[2];
    const double *Sx = in[3];
    const double *Su = in[4];

    double *jac_x = out[0];
    double *jac_u = out[1];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[5];
    double *casadi_res[2];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;
    casadi_arg[3] = Sx;
    casadi_arg[4] = Su;

    casadi_res[0] = jac_x;
    casadi_res[1] = jac_u;

    vdeFun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void impl_f_Fun(double **in, double **out){
    const double *x = in[0];
    const double *u = in[1];
    const double *p = in[2];
    const double *xdot = in[3];

    double *res = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;
    casadi_arg[3] = xdot;

    casadi_res[0] = res;

    impl_f_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void impl_jac_x_Fun(double **in, double **out){
    const double *x = in[0];
    const double *u = in[1];
    const double *p = in[2];
    const double *xdot = in[3];

    double *jac_x = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;
    casadi_arg[3] = xdot;

    casadi_res[0] = jac_x;

    impl_jac_x_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void impl_jac_u_Fun(double **in, double **out){
    const double *x = in[0];
    const double *u = in[1];
    const double *p = in[2];
    const double *xdot = in[3];

    double *jac_u = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;
    casadi_arg[3] = xdot;

    casadi_res[0] = jac_u;

    impl_jac_u_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void impl_jac_xdot_Fun(double **in, double **out){
    const double *x = in[0];
    const double *u = in[1];
    const double *p = in[2];
    const double *xdot = in[3];

    double *jac_xdot = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;
    casadi_arg[3] = xdot;

    casadi_res[0] = jac_xdot;

    impl_jac_xdot_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void F_Fun(double **in, double **out) {

    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];

    double *x_out = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;

    casadi_res[0] = x_out;

    F(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void D_Fun(double **in, double **out) {

    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];

    double *A = out[0];
    double *B = out[1];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[2];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;

    casadi_res[0] = A;
    casadi_res[1] = B;

    D(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void Hi_Fun(double **in, double **out) {
    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];
    const double *refi = in[3];
    const double *Qi = in[4];

    double *Hx = out[0];
    double *Hu = out[1];
    double *Hxu = out[2];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[5];
    double *casadi_res[3];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;
    casadi_arg[3] = refi;
    casadi_arg[4] = Qi;

    casadi_res[0] = Hx;
    casadi_res[1] = Hu;
    casadi_res[2] = Hxu;

    Hi_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void HN_Fun(double **in, double **out) {
    const double *xN = in[0];
    const double *paraN = in[1] ;
    const double *refN = in[2];
    const double *QN = in[3];

    double *HN = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[1];

    casadi_arg[0] = xN;
    casadi_arg[1] = paraN;
    casadi_arg[2] = refN;
    casadi_arg[3] = QN;

    casadi_res[0] = HN;

    HN_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void gi_Fun(double **in, double **out) {

    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];
    const double *refi = in[3];
    const double *Qi = in[4];

    double *gx = out[0];
    double *gu = out[1];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[5];
    double *casadi_res[2];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;
    casadi_arg[3] = refi;
    casadi_arg[4] = Qi;

    casadi_res[0] = gx;
    casadi_res[1] = gu;

    gi_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void path_con_Fun(double **in, double **out) {

    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];

    double *ci = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;

    casadi_res[0] = ci;

    path_con_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void Ci_Fun(double **in, double **out) {
    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];

    double *Cx = out[0];
    double *Cu = out[1];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[2];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;

    casadi_res[0] = Cx;
    casadi_res[1] = Cu;

    Ci_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void gN_Fun(double **in, double **out) {

    const double *xN = in[0];
    const double *paraN = in[1] ;
    const double *refN = in[2];
    const double *QN = in[3];

    double *gxN = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[1];

    casadi_arg[0] = xN;
    casadi_arg[1] = paraN;
    casadi_arg[2] = refN;
    casadi_arg[3] = QN;

    casadi_res[0] = gxN;

    gN_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void path_con_N_Fun(double **in, double **out) {

    const double *xN = in[0];
    const double *paraN = in[1];

    double *cN = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[2];
    double *casadi_res[1];

    casadi_arg[0] = xN;
    casadi_arg[1] = paraN;

    casadi_res[0] = cN;

    path_con_N_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void CN_Fun(double **in, double **out) {

    const double *xN = in[0];
    const double *paraN = in[1];

    double *CxN = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[2];
    double *casadi_res[1];

    casadi_arg[0] = xN;
    casadi_arg[1] = paraN;

    casadi_res[0] = CxN;

    CN_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void adj_Fun(double **in, double **out) {

    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];
    const double *refi = in[3];
    const double *Qi = in[4];
    const double *lambdai = in[5];
    const double *mui = in[6];

    double *dobj = out[0];
    double *adj_dG = out[1];
    double *adj_dB = out[2];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[7];
    double *casadi_res[3];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;
    casadi_arg[3] = refi;
    casadi_arg[4] = Qi;
    casadi_arg[5] = lambdai;
    casadi_arg[6] = mui;

    casadi_res[0] = dobj;
    casadi_res[1] = adj_dG;
    casadi_res[2] = adj_dB;

    adj_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void adjN_Fun(double **in, double **out) {

    const double *xN = in[0];
    const double *paraN = in[1];
    const double *refN = in[2];
    const double *QN = in[3];
    const double *muN = in[4];

    double *dobj = out[0];
    double *adj_dB = out[1];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[5];
    double *casadi_res[2];

    casadi_arg[0] = xN;
    casadi_arg[1] = paraN;
    casadi_arg[2] = refN;
    casadi_arg[3] = QN;
    casadi_arg[4] = muN;

    casadi_res[0] = dobj;
    casadi_res[1] = adj_dB;

    adjN_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void obji_Fun(double **in, double **out){
    const double *xi = in[0];
    const double *ui = in[1];
    const double *parai = in[2];
    const double *refi = in[3];
    const double *Qi = in[4];

    double *obj = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[5];
    double *casadi_res[1];

    casadi_arg[0] = xi;
    casadi_arg[1] = ui;
    casadi_arg[2] = parai;
    casadi_arg[3] = refi;
    casadi_arg[4] = Qi;

    casadi_res[0] = obj;

    obji_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void objN_Fun(double **in, double **out){
    const double *xN = in[0];
    const double *parai = in[1];
    const double *refN = in[2];
    const double *QN = in[3];

    double *objN = out[0];

    int casadi_mem = 0;
    long long int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[1];

    casadi_arg[0] = xN;
    casadi_arg[1] = parai;
    casadi_arg[2] = refN;
    casadi_arg[3] = QN;

    casadi_res[0] = objN;

    objN_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}
