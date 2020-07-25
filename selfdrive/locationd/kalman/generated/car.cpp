
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4011336923141012289) {
   out_4011336923141012289[0] = delta_x[0] + nom_x[0];
   out_4011336923141012289[1] = delta_x[1] + nom_x[1];
   out_4011336923141012289[2] = delta_x[2] + nom_x[2];
   out_4011336923141012289[3] = delta_x[3] + nom_x[3];
   out_4011336923141012289[4] = delta_x[4] + nom_x[4];
   out_4011336923141012289[5] = delta_x[5] + nom_x[5];
   out_4011336923141012289[6] = delta_x[6] + nom_x[6];
   out_4011336923141012289[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8823850480828624609) {
   out_8823850480828624609[0] = -nom_x[0] + true_x[0];
   out_8823850480828624609[1] = -nom_x[1] + true_x[1];
   out_8823850480828624609[2] = -nom_x[2] + true_x[2];
   out_8823850480828624609[3] = -nom_x[3] + true_x[3];
   out_8823850480828624609[4] = -nom_x[4] + true_x[4];
   out_8823850480828624609[5] = -nom_x[5] + true_x[5];
   out_8823850480828624609[6] = -nom_x[6] + true_x[6];
   out_8823850480828624609[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_987266334863380994) {
   out_987266334863380994[0] = 1.0;
   out_987266334863380994[1] = 0.0;
   out_987266334863380994[2] = 0.0;
   out_987266334863380994[3] = 0.0;
   out_987266334863380994[4] = 0.0;
   out_987266334863380994[5] = 0.0;
   out_987266334863380994[6] = 0.0;
   out_987266334863380994[7] = 0.0;
   out_987266334863380994[8] = 0.0;
   out_987266334863380994[9] = 1.0;
   out_987266334863380994[10] = 0.0;
   out_987266334863380994[11] = 0.0;
   out_987266334863380994[12] = 0.0;
   out_987266334863380994[13] = 0.0;
   out_987266334863380994[14] = 0.0;
   out_987266334863380994[15] = 0.0;
   out_987266334863380994[16] = 0.0;
   out_987266334863380994[17] = 0.0;
   out_987266334863380994[18] = 1.0;
   out_987266334863380994[19] = 0.0;
   out_987266334863380994[20] = 0.0;
   out_987266334863380994[21] = 0.0;
   out_987266334863380994[22] = 0.0;
   out_987266334863380994[23] = 0.0;
   out_987266334863380994[24] = 0.0;
   out_987266334863380994[25] = 0.0;
   out_987266334863380994[26] = 0.0;
   out_987266334863380994[27] = 1.0;
   out_987266334863380994[28] = 0.0;
   out_987266334863380994[29] = 0.0;
   out_987266334863380994[30] = 0.0;
   out_987266334863380994[31] = 0.0;
   out_987266334863380994[32] = 0.0;
   out_987266334863380994[33] = 0.0;
   out_987266334863380994[34] = 0.0;
   out_987266334863380994[35] = 0.0;
   out_987266334863380994[36] = 1.0;
   out_987266334863380994[37] = 0.0;
   out_987266334863380994[38] = 0.0;
   out_987266334863380994[39] = 0.0;
   out_987266334863380994[40] = 0.0;
   out_987266334863380994[41] = 0.0;
   out_987266334863380994[42] = 0.0;
   out_987266334863380994[43] = 0.0;
   out_987266334863380994[44] = 0.0;
   out_987266334863380994[45] = 1.0;
   out_987266334863380994[46] = 0.0;
   out_987266334863380994[47] = 0.0;
   out_987266334863380994[48] = 0.0;
   out_987266334863380994[49] = 0.0;
   out_987266334863380994[50] = 0.0;
   out_987266334863380994[51] = 0.0;
   out_987266334863380994[52] = 0.0;
   out_987266334863380994[53] = 0.0;
   out_987266334863380994[54] = 1.0;
   out_987266334863380994[55] = 0.0;
   out_987266334863380994[56] = 0.0;
   out_987266334863380994[57] = 0.0;
   out_987266334863380994[58] = 0.0;
   out_987266334863380994[59] = 0.0;
   out_987266334863380994[60] = 0.0;
   out_987266334863380994[61] = 0.0;
   out_987266334863380994[62] = 0.0;
   out_987266334863380994[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5568542721591178760) {
   out_5568542721591178760[0] = state[0];
   out_5568542721591178760[1] = state[1];
   out_5568542721591178760[2] = state[2];
   out_5568542721591178760[3] = state[3];
   out_5568542721591178760[4] = state[4];
   out_5568542721591178760[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5568542721591178760[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5568542721591178760[7] = state[7];
}
void F_fun(double *state, double dt, double *out_1865630451206577316) {
   out_1865630451206577316[0] = 1;
   out_1865630451206577316[1] = 0;
   out_1865630451206577316[2] = 0;
   out_1865630451206577316[3] = 0;
   out_1865630451206577316[4] = 0;
   out_1865630451206577316[5] = 0;
   out_1865630451206577316[6] = 0;
   out_1865630451206577316[7] = 0;
   out_1865630451206577316[8] = 0;
   out_1865630451206577316[9] = 1;
   out_1865630451206577316[10] = 0;
   out_1865630451206577316[11] = 0;
   out_1865630451206577316[12] = 0;
   out_1865630451206577316[13] = 0;
   out_1865630451206577316[14] = 0;
   out_1865630451206577316[15] = 0;
   out_1865630451206577316[16] = 0;
   out_1865630451206577316[17] = 0;
   out_1865630451206577316[18] = 1;
   out_1865630451206577316[19] = 0;
   out_1865630451206577316[20] = 0;
   out_1865630451206577316[21] = 0;
   out_1865630451206577316[22] = 0;
   out_1865630451206577316[23] = 0;
   out_1865630451206577316[24] = 0;
   out_1865630451206577316[25] = 0;
   out_1865630451206577316[26] = 0;
   out_1865630451206577316[27] = 1;
   out_1865630451206577316[28] = 0;
   out_1865630451206577316[29] = 0;
   out_1865630451206577316[30] = 0;
   out_1865630451206577316[31] = 0;
   out_1865630451206577316[32] = 0;
   out_1865630451206577316[33] = 0;
   out_1865630451206577316[34] = 0;
   out_1865630451206577316[35] = 0;
   out_1865630451206577316[36] = 1;
   out_1865630451206577316[37] = 0;
   out_1865630451206577316[38] = 0;
   out_1865630451206577316[39] = 0;
   out_1865630451206577316[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1865630451206577316[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1865630451206577316[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1865630451206577316[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1865630451206577316[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1865630451206577316[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1865630451206577316[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1865630451206577316[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1865630451206577316[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1865630451206577316[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1865630451206577316[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1865630451206577316[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1865630451206577316[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1865630451206577316[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1865630451206577316[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1865630451206577316[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1865630451206577316[56] = 0;
   out_1865630451206577316[57] = 0;
   out_1865630451206577316[58] = 0;
   out_1865630451206577316[59] = 0;
   out_1865630451206577316[60] = 0;
   out_1865630451206577316[61] = 0;
   out_1865630451206577316[62] = 0;
   out_1865630451206577316[63] = 1;
}
void h_25(double *state, double *unused, double *out_7910287640668498307) {
   out_7910287640668498307[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4423232291322955318) {
   out_4423232291322955318[0] = 0;
   out_4423232291322955318[1] = 0;
   out_4423232291322955318[2] = 0;
   out_4423232291322955318[3] = 0;
   out_4423232291322955318[4] = 0;
   out_4423232291322955318[5] = 0;
   out_4423232291322955318[6] = 1;
   out_4423232291322955318[7] = 0;
}
void h_24(double *state, double *unused, double *out_314578370085193696) {
   out_314578370085193696[0] = state[4];
   out_314578370085193696[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1447805804439264326) {
   out_1447805804439264326[0] = 0;
   out_1447805804439264326[1] = 0;
   out_1447805804439264326[2] = 0;
   out_1447805804439264326[3] = 0;
   out_1447805804439264326[4] = 1;
   out_1447805804439264326[5] = 0;
   out_1447805804439264326[6] = 0;
   out_1447805804439264326[7] = 0;
   out_1447805804439264326[8] = 0;
   out_1447805804439264326[9] = 0;
   out_1447805804439264326[10] = 0;
   out_1447805804439264326[11] = 0;
   out_1447805804439264326[12] = 0;
   out_1447805804439264326[13] = 1;
   out_1447805804439264326[14] = 0;
   out_1447805804439264326[15] = 0;
}
void h_30(double *state, double *unused, double *out_5620735426563764127) {
   out_5620735426563764127[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3100318507738811744) {
   out_3100318507738811744[0] = 0;
   out_3100318507738811744[1] = 0;
   out_3100318507738811744[2] = 0;
   out_3100318507738811744[3] = 0;
   out_3100318507738811744[4] = 1;
   out_3100318507738811744[5] = 0;
   out_3100318507738811744[6] = 0;
   out_3100318507738811744[7] = 0;
}
void h_26(double *state, double *unused, double *out_2329839985467055884) {
   out_2329839985467055884[0] = state[7];
}
void H_26(double *state, double *unused, double *out_119992176909865698) {
   out_119992176909865698[0] = 0;
   out_119992176909865698[1] = 0;
   out_119992176909865698[2] = 0;
   out_119992176909865698[3] = 0;
   out_119992176909865698[4] = 0;
   out_119992176909865698[5] = 0;
   out_119992176909865698[6] = 0;
   out_119992176909865698[7] = 1;
}
void h_27(double *state, double *unused, double *out_3320923911940229360) {
   out_3320923911940229360[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2460380156234230334) {
   out_2460380156234230334[0] = 0;
   out_2460380156234230334[1] = 0;
   out_2460380156234230334[2] = 0;
   out_2460380156234230334[3] = 1;
   out_2460380156234230334[4] = 0;
   out_2460380156234230334[5] = 0;
   out_2460380156234230334[6] = 0;
   out_2460380156234230334[7] = 0;
}
void h_29(double *state, double *unused, double *out_8320824043877247602) {
   out_8320824043877247602[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8565204038721919226) {
   out_8565204038721919226[0] = 0;
   out_8565204038721919226[1] = 1;
   out_8565204038721919226[2] = 0;
   out_8565204038721919226[3] = 0;
   out_8565204038721919226[4] = 0;
   out_8565204038721919226[5] = 0;
   out_8565204038721919226[6] = 0;
   out_8565204038721919226[7] = 0;
}
void h_28(double *state, double *unused, double *out_5425862456618261424) {
   out_5425862456618261424[0] = state[5];
   out_5425862456618261424[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7836969225713155748) {
   out_7836969225713155748[0] = 0;
   out_7836969225713155748[1] = 0;
   out_7836969225713155748[2] = 0;
   out_7836969225713155748[3] = 0;
   out_7836969225713155748[4] = 0;
   out_7836969225713155748[5] = 1;
   out_7836969225713155748[6] = 0;
   out_7836969225713155748[7] = 0;
   out_7836969225713155748[8] = 0;
   out_7836969225713155748[9] = 0;
   out_7836969225713155748[10] = 0;
   out_7836969225713155748[11] = 0;
   out_7836969225713155748[12] = 0;
   out_7836969225713155748[13] = 0;
   out_7836969225713155748[14] = 1;
   out_7836969225713155748[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
