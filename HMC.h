
#include "mix_lib.h"

extern Generator g;
extern int Kmax;
extern int Curve_num;
extern int Nm;
extern double w_eps,v0_eps,sigmav2_eps;
extern ofstream testfile;

double K(const pq_point &theta, int k, const double lambda);
double U(const curve Data[], const pq_point &theta, int k, const vec & z);
double logP_partial_w_percurve_cpp(const curve &Dat, const pq_point &theta, int k);
vec U_partial_w_cpp(const curve Data[], const pq_point &theta, const vec & z);
double logP_partial_v_percurve_cpp(const curve &Dat, const pq_point &theta, int k);
vec U_partial_v_cpp(const curve Data[], const pq_point &theta, const vec & z);
double logP_partial_sigma2_percurve_cpp(const curve &Dat, const pq_point &theta, int k);
vec U_partial_sigma2_cpp(const curve Data[],const pq_point &theta, const vec & z);
void Gibbs_Sampling_pi( pq_point &theta, const vec & z);
void Gibbs_Sampling_z(const curve Data[], const pq_point & theta, vec & z);
