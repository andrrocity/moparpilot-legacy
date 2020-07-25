/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7299192101197356359);
void inv_err_fun(double *nom_x, double *true_x, double *out_9089171002813238635);
void H_mod_fun(double *state, double *out_258406204140633298);
void f_fun(double *state, double dt, double *out_8495819442745364591);
void F_fun(double *state, double dt, double *out_4727433309476161094);
void h_3(double *state, double *unused, double *out_7649419506753245627);
void H_3(double *state, double *unused, double *out_4752749176699099027);
void h_4(double *state, double *unused, double *out_2633841578518801850);
void H_4(double *state, double *unused, double *out_1700867757280659926);
void h_9(double *state, double *unused, double *out_6170250320230314691);
void H_9(double *state, double *unused, double *out_1395333267518808676);
void h_10(double *state, double *unused, double *out_6808254535763801354);
void H_10(double *state, double *unused, double *out_5000233975320603562);
void h_12(double *state, double *unused, double *out_4468212991041160513);
void H_12(double *state, double *unused, double *out_2534005910367908180);
void h_13(double *state, double *unused, double *out_7887645412769034432);
void H_13(double *state, double *unused, double *out_2911383651817869696);
void h_14(double *state, double *unused, double *out_6170250320230314691);
void H_14(double *state, double *unused, double *out_1395333267518808676);
void h_19(double *state, double *unused, double *out_7785984444803071347);
void H_19(double *state, double *unused, double *out_7691230213154749720);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);