/********************************************************************************************************************/
/*mix_lib.h: Header file for mix_lib, rj_mix.   */
/********************************************************************************************************************/

#define MaxM 2000
#define MaxL 102
#define MaxK 30
// #define RAND_MAX 0x7fff
#define NSample 100000

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

#include "alea.h"
#include <armadillo>
#include <iostream>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <omp.h>
using namespace std;
using namespace arma;

/********************************************************************************************************************/
/*Useful functions aliases and definitions */
/********************************************************************************************************************/
#define NEG_REAL_MIN (-1.79e+308)
#define eps 1e-7

using namespace arma;
using namespace std;
/********************************************************************************************************************/
/*Data types */
/********************************************************************************************************************/
struct pq_point
{
        vec pi;
        vec w;
        vec v;
        vec sigma2;

        explicit pq_point(int K) : pi(K), w(K), v(K), sigma2(K)
        {
        }
        pq_point()
        {
        }
        pq_point(const pq_point &theta) : pi(theta.pi.size()), w(theta.w.size()), v(theta.v.size()), sigma2(theta.sigma2.size())
        {
                pi = theta.pi;
                w = theta.w;
                v = theta.v;
                sigma2 = theta.sigma2;
        }

        pq_point &operator=(const pq_point &theta)
        {
                if (this == &theta)
                        return *this;

                pi = theta.pi;
                w = theta.w;
                v = theta.v;
                sigma2 = theta.sigma2;

                return *this;
        }
        pq_point operator+(const pq_point &theta) const
        {
                pq_point sum;
                sum.w = w + theta.w;
                sum.v = v + theta.v;
                sum.sigma2 = sigma2 + theta.sigma2;

                return sum;
        }

        pq_point operator%(const pq_point &theta) const
        {
                pq_point multi;
                multi.w = w % theta.w;
                multi.v = v % theta.v;
                multi.sigma2 = sigma2 % theta.sigma2;

                return multi;
        }

        double accu() const
        {
                return sum(w) + sum(v) + sum(sigma2);
        }

        void zeros()
        {
                w.zeros();
                v.zeros();
                sigma2.zeros();
        }

        void initialize(Generator &g)
        {
                int K = w.size();

                for (int k = 0; k < K; k++)
                {
                        w(k) = 0.5 / Gam(0.5, g);
                        v(k) = 0.5 / Gam(0.5, g);
                        sigma2(k) = 0.01 / Gam(10, g);
                }
        }

        void positive_reflect()
        {
                int K = w.size();
                for (int k = 0; k < K; k++)
                {
                        if (w(k) < 0)
                                w(k) = -w(k);
                        if (v(k) < 0)
                                v(k) = -v(k);
                        if (sigma2(k) < 0)
                                sigma2(k) = -sigma2(k);
                }
        }

        void draw_from_prior(Generator &g)
        {
                // double ran1,ran2;
                int K = w.size();
                // cout << K << endl;
                for (int k = 0; k < K; k++)
                {
                        double tmp = kiss(g);
                        // cout<< tmp<<" ";
                        pi(k) = tmp;
                        w(k) = 0.5 / Gam(0.5, g);
                        v(k) = 0.5 / Gam(0.5, g);
                        sigma2(k) = 0.01 / Gam(10, g);
                }
                // cout<<endl;
                // pi.print("before normalise pi:");
                pi = pi/sum(pi);
                // pi.print("after normalise pi:");

                // *this.print("theta:");
        }

        int locate_seq(double q)
        {
                for (int k = 0; k < w.size(); k++)
                        if (w(k) > q)
                                return k;

                return w.size();
        }

        void insertPre_seq(int k, pq_point Elem_mix)
        {
                int K = w.size();

                if (k < 0)
                { //?k>K?
                        pi.insert_rows(K, Elem_mix.pi);
                        w.insert_rows(K, Elem_mix.w);
                        v.insert_rows(K, Elem_mix.v);
                        sigma2.insert_rows(K, Elem_mix.sigma2);
                }

                else
                {
                        pi.insert_rows(k, Elem_mix.pi);
                        w.insert_rows(k, Elem_mix.w);
                        v.insert_rows(k, Elem_mix.v);
                        sigma2.insert_rows(k, Elem_mix.sigma2);
                }
        }
        void deleteP_seq(int k)
        {
                pi.shed_row(k);
                w.shed_row(k);
                v.shed_row(k);
                sigma2.shed_row(k);
        }
        void deleteP_seq(uvec &empty_component)
        {
                int k1 = *(empty_component.begin_row(1)), k2 = *(empty_component.end_row(1));
                pi.shed_rows(k1, k2);
                w.shed_rows(k1, k2);
                v.shed_rows(k1, k2);
                sigma2.shed_rows(k1, k2);
        }
        void print() const
        {
                cout << "K: " << pi.size() << endl;
                pi.print("pi:");
                w.print("w:");
                v.print("v:");
                sigma2.print("sigma2:");
        }

        void print(string s) const
        {
                cout << s << endl;
                cout << "K: " << pi.size() << endl;
                pi.print("pi:");
                w.print("w:");
                v.print("v:");
                sigma2.print("sigma2:");
        }
        void write_file(ofstream &myfile)
        {
                myfile << pi << w << v << sigma2;
        }
        void write_file(string s, ofstream &myfile)
        {
                myfile << s << endl;
                myfile << pi << w << v << sigma2;
        }
};

struct curve
{

        vec X;
        vec Y;
        explicit curve(int Nm) : X(Nm), Y(Nm)
        {
        }
        curve()
        {
        }
        curve(const curve &Dat) : X(Dat.X.size()), Y(Dat.Y.size())
        {
                X = Dat.X;
                Y = Dat.Y;
        }
        void write_file(int i, ofstream &myfile)
        {
                myfile << i << "th curve data:" << endl;
                myfile << "X is:" << X << "Y is:" << Y;
        }
};

/********************************************************************************************************************/
/*global variables */
/********************************************************************************************************************/
/* Seeds for the kiss generator (see alea) */

/*Model and sampler parameters(defined in main program)  */
extern Generator g;
extern int Kmax;
extern int Curve_num;
extern int Nm;
extern double w_eps, v0_eps, sigmav2_eps;
extern ofstream testfile;

/********************************************************************************************************************/
/*Routines. */
/********************************************************************************************************************/
void xixj(mat & K, const vec &X);
double log_likelihood_micro2(const curve &Dat, const pq_point &theta, int k);
double log_likelihood2(const curve Data[], const pq_point &theta);
double logSumExp(const vec& x);
void draw_initial_model(const curve Data[], pq_point &theta, double *logl);
double prop_split(const curve Data[], pq_point &theta, int k, int *k1, int *k2);
double prop_merge(const curve Data[], pq_point &theta_old, pq_point &theta_new, int *k, int k1, int k2);
double compute_log_split_ratio(pq_point &theta, pq_point &theta_split, int k, int k_split1, int k_split2);
double calc_secondary_moment(const curve &Dat, const pq_point &theta, const pq_point &theta_new, int k, int k1, int k2);
double compute_log_prior_ratio(const pq_point &theta, const pq_point &theta_old);
void zeros(int c[], int K);
