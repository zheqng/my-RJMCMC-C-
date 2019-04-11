
#include "HMC.h"

double K(const pq_point &theta, int k, const double lambda)
{
        // cout<<"calc K"<<endl;
        vec kinetic_energy = zeros<vec>(1);
        kinetic_energy += pow(theta.w(k), 2.0) / 2.0 / lambda + pow(theta.v(k), 2.0) / 2.0 / lambda + pow(theta.sigma2(k), 2.0) / 2.0 / lambda;
        return kinetic_energy(1);
}

double U(const curve Data[], const pq_point &theta, int k, const vec &z)
{
        // cout<<" calc U"<<endl;
        double u = 0.0, l_pk = 0.0;
        for (int m = 0; m < Curve_num; m++)
        {
                if (z(m) == k)
                {
                        l_pk += log_likelihood_micro2(Data[m], theta, k);
                }
        }
        boost::math::inverse_gamma_distribution<> IG_dist_w(0.5, 0.5);
        boost::math::inverse_gamma_distribution<> IG_dist_v(0.5, 0.5);
        boost::math::inverse_gamma_distribution<> IG_dist_s(10, 0.01);
        // theta.print("theta:");
        double lp_w = log(pdf(IG_dist_w, theta.w(k)));
        double lp_v = log(pdf(IG_dist_v, theta.v(k)));
        double lp_sigma2 = log(pdf(IG_dist_s, theta.sigma2(k)));

        u = -l_pk - lp_w - lp_v - lp_sigma2;
        return u;
}

/*U0 partial w*/
double logP_partial_w_percurve_cpp(const curve &Dat, const pq_point &theta, int k)
{
        int r = Dat.X.size();
        mat cov(r, r), covinv(r, r), cov_partial_w(r, r), squareX(r, r);
        vec tmp(1);
        vec covinv_Y(r);

        // vec colones = ones<vec>(r);
        // rowvec rowones = ones<rowvec>(r);
        mat EYE = eye<mat>(r, r);
        xixj(squareX,Dat.X);
        // square(Dat.X * rowones - colones * trans(Dat.X));
        cov = theta.v(k) * exp(-squareX * theta.w(k) / 2) + EYE * theta.sigma2(k);
        // covinv = inv(cov);
        cov_partial_w = -0.5 * theta.v(k) * (squareX % exp(-squareX * theta.w(k) / 2.0));
        // tmp = -0.5 * trace(covinv * cov_partial_w);
        tmp = -0.5 * trace(solve(cov, cov_partial_w));
        covinv_Y = solve(cov, Dat.Y);
        // tmp = tmp + 0.5 * trans(Dat.Y) * covinv * cov_partial_w * covinv * Dat.Y;
        tmp = tmp + 0.5 * trans(covinv_Y) * cov_partial_w * covinv_Y;
        return tmp(0);
}

vec U_partial_w_cpp(const curve Data[], const pq_point &theta, const vec &z)
{
        int K = theta.w.size();
        vec u_partial_w = zeros<vec>(K);
        for (int k = 0; k < K; k++)
        {
                for (int m = 0; m < Curve_num; m++)
                {
                        if (z(m) == k)
                        {

                                u_partial_w(k) = u_partial_w(k) - logP_partial_w_percurve_cpp(Data[m], theta, k);
                        }
                }
        }
        u_partial_w = u_partial_w + (0.5 + 1) / theta.w - 0.5 / square(theta.w);

        return u_partial_w;
}
/*U0 partial v0*/
double logP_partial_v_percurve_cpp(const curve &Dat, const pq_point &theta, int k)
{
        int r = Dat.X.size();
        mat cov(r, r), covinv(r, r), cov_partial_v(r, r), squareX(r, r);
        vec tmp(1);
        vec covinv_Y(r);

        // vec colones = ones<vec>(r);
        // rowvec rowones = ones<rowvec>(r);
        mat EYE = eye<mat>(r, r);
        xixj(squareX,Dat.X);
        // square(Dat.X * rowones - colones * trans(Dat.X));

        cov = theta.v(k) * exp(-squareX * theta.w(k) / 2.0) + EYE * theta.sigma2(k);
        // covinv = inv(cov);
        cov_partial_v = exp(-squareX * theta.w(k) / 2.0);
        covinv_Y = solve(cov, Dat.Y);
        tmp = -0.5 * trace(solve(cov, cov_partial_v)) + 0.5 * trans(covinv_Y) * cov_partial_v * covinv_Y;

        return tmp(0);
}

vec U_partial_v_cpp(const curve Data[], const pq_point &theta, const vec &z)
{
        int K = theta.w.size();
        vec u_partial_v = zeros<vec>(K);
        for (int k = 0; k < K; k++)
        {
                for (int m = 0; m < Curve_num; m++)
                {
                        if (z(m) == k)
                        {
                                u_partial_v(k) = u_partial_v(k) - logP_partial_v_percurve_cpp(Data[m], theta, k);
                        }
                }
        }
        u_partial_v += (0.5 + 1) / theta.v - 0.5 / square(theta.v);
        return u_partial_v;
}

/*U0 partial sigmav2*/
double logP_partial_sigma2_percurve_cpp(const curve &Dat, const pq_point &theta, int k)
{
        int r = Dat.X.size();
        mat cov(r, r), covinv(r, r), cov_partial_sigma2(r, r), squareX(r, r);
        vec tmp(1);
        vec covinv_Y(r);
        // caculate cov
        // vec colones = ones<vec>(r);
        // rowvec rowones = ones<rowvec>(r);
        mat EYE = eye<mat>(r, r);
        xixj(squareX,Dat.X);
        // square(Dat.X * rowones - colones * trans(Dat.X));

        cov = theta.v(k) * exp(-squareX * theta.w(k) / 2.0) + EYE * theta.sigma2(k);
        // covinv = inv(cov);
        cov_partial_sigma2 = EYE;
        covinv_Y = solve(cov,Dat.Y);
        tmp = -0.5 * trace( solve(cov, cov_partial_sigma2)) + 0.5 * trans(covinv_Y) *  cov_partial_sigma2 * covinv_Y;

        return tmp(0);
}

vec U_partial_sigma2_cpp(const curve Data[], const pq_point &theta, const vec &z)
{
        int K = theta.w.size();
        vec u_partial_sigma2 = zeros<vec>(K);
        for (int m = 0; m < Curve_num; m++)
        {
                for (int k = 0; k < K; k++)
                {
                        if (z(m) == k)
                        {
                                u_partial_sigma2(k) = u_partial_sigma2(k) - logP_partial_sigma2_percurve_cpp(Data[m], theta, k);
                        }
                }
        }
        u_partial_sigma2 += (0.5 + 1) / theta.sigma2 - 0.5 / square(theta.sigma2);
        return u_partial_sigma2;
}

/*HMC */

/*______________________________Gibbs Sampling____________________________*/

void Gibbs_Sampling_pi(pq_point &theta, const vec &z)
{
        int c[MaxK];
        double sum = 0.0;
        int K = theta.w.size();

        zeros(c, K);
        for (int n = 0; n < Curve_num; n++)
        {
                int k = z(n);
                c[k] = c[k] + 1;
        }

        for (int k = 0; k < theta.w.size(); k++)
        {
                theta.pi(k) = Gam(1 + c[k], g);
        }

        theta.pi = normalise(theta.pi, 1);
}

void Gibbs_Sampling_z(const curve Data[], const pq_point &theta, vec &z)
{
        int K = theta.w.size();
        vec logp_k(K);
        z = zeros<vec>(Curve_num);
        for (int m = 0; m < Curve_num; m++)
        {
                for (int k = 0; k < K; k++)
                {
                        logp_k(k) =log_likelihood_micro2(Data[m], theta, k);
                }
                /*_________________________hard cut___________________________*/
                logp_k = logp_k +  log(theta.pi);
                z(m) = logp_k.index_max();
        }
}
