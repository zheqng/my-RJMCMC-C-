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
#include <Eigen/Dense>
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <omp.h>
// using namespace std;
using namespace Eigen;

/********************************************************************************************************************/
/*Useful functions aliases and definitions */
/********************************************************************************************************************/
#define NEG_REAL_MIN (-1.79e+308)
// #define eps 1e-7

// using namespace arma;
// using namespace std;
extern Generator g;
extern int Kmax;
extern int Curve_num, THREAD_NUM;
extern int Nm;
extern double w_eps, v0_eps, sigmav2_eps;
extern ofstream testfile;

/********************************************************************************************************************/
/*Data types */
/********************************************************************************************************************/
struct pq_point
{
        VectorXd pi;
        VectorXd w;
        VectorXd v;
        VectorXd sigma2;

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
        double operator%(const pq_point &theta) const
        {
                double multi;
                // cout<<"operator%****************"<<endl;
                // cout<<"w:"<<w<<endl;
                // cout<<"theta.w:"<<theta.w<<endl;
                // cout<<" w.dot(theta.w):"<< w.dot(theta.w)<<"  v.dot(theta.v):"<<v.dot(theta.v)<<"  sigma2.dot(theta.sigma2):"<<sigma2.dot(theta.sigma2)<<endl;

                multi = w.dot(theta.w) +v.dot(theta.v) +sigma2.dot(theta.sigma2);
                // multi.v = v.dot(theta.v);
                // multi.sigma2 = sigma2.dot(theta.sigma2);

                return multi;
        }
        double sum() const
        {
                return w.sum() + v.sum() + sigma2.sum();
        }



        void zeros()
        {
                int K = w.size();
                w  =  VectorXd::Zero(K);
                v =  VectorXd::Zero(K);
                sigma2 =  VectorXd::Zero(K);
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
                pi = pi/pi.sum();
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
                // cout<<"k="<<k<<" K="<<K<<endl;
                // cout<<v<<endl;
                VectorXd tmp = pi; pi.resize(K+1); pi.head(K) = tmp;
                tmp = w; w.resize(K+1); w.head(K) = tmp;
                tmp = v; v.resize(K+1); v.head(K) = tmp;
                tmp = sigma2; sigma2.resize(K+1); sigma2.head(K) = tmp;
                if (k >=K)
                {   //?k>K?
                    //         pi.insert_rows(K, Elem_mix.pi);
                    //         w.insert_rows(K, Elem_mix.w);
                    //         v.insert_rows(K, Elem_mix.v);
                    //         sigma2.insert_rows(K, Elem_mix.sigma2);
                        pi(k) = Elem_mix.pi(0);
                        w(k) = Elem_mix.w(0);
                        v(k) = Elem_mix.v(0);
                        sigma2(k) = Elem_mix.sigma2(0);
                        // cout<<v<<endl;
                }

                else
                {
                        pi.tail(K - k) = pi.segment(k, K - 1); pi(k) =Elem_mix.pi(0);
                        w.tail(K - k) = w.segment(k, K - 1); w(k) =Elem_mix.w(0);
                        v.tail(K - k) = v.segment(k, K - 1); v(k) =Elem_mix.v(0);
                        sigma2.tail(K - k) = sigma2.segment(k, K - 1); sigma2(k) =Elem_mix.sigma2(0);
                }
                // cout<<"after insert_seq, the size is accumulated by:"<<pi.size() - K<<endl;

        }
        void deleteP_seq(int k)
        {
                int K = pi.size();
                VectorXd tmp(K-1); tmp.head(k) = pi.head(k);
                if(K==k) {}
                else{tmp.tail(K-k-1)= pi.tail(K-k-1);}
                pi.resize(K-1); pi=tmp;

                tmp.head(k) = w.head(k);
                if(K==k) {}
                else{tmp.tail(K-k-1)= w.tail(K-k-1);}
                w.resize(K-1); w=tmp;

                // cout<<"v:"<<v<<endl;
                tmp.head(k) = v.head(k);
                if(K==k) {}
                else{tmp.tail(K-k-1)= v.tail(K-k-1);}
                v.resize(K-1); v=tmp;
                // cout<<"v:"<<v<<endl;

                tmp.head(k) = sigma2.head(k);
                if(K==k) {}
                else{tmp.tail(K-k-1)= sigma2.tail(K-k-1);}
                sigma2.resize(K-1); sigma2=tmp;

        }

        void print() const
        {
                cout << "K: " << pi.size() << endl;
                cout<<"pi:\n"<<pi<<endl;;
                cout<<"w:\n"<<w<<endl;
                cout<<"v:\n"<<v<<endl;
                cout<<"sigma2:\n"<<sigma2<<endl;
        }

        void print(string s) const
        {
                cout << s << endl;
                cout << "K: " << pi.size() << endl;
                cout<<"pi:\n"<<pi<<endl;;
                cout<<"w:\n"<<w<<endl;
                cout<<"v:\n"<<v<<endl;
                cout<<"sigma2:\n"<<sigma2<<endl;
        }
        void write_file(ofstream &myfile)
        {
                myfile << pi <<endl;
                myfile<< w << endl;
                myfile << v <<endl;
                myfile << sigma2<<endl;
        }
        void write_file(string s, ofstream &myfile)
        {
                myfile << s << endl;
                myfile << pi <<endl;
                myfile<< w << endl;
                myfile << v <<endl;
                myfile << sigma2<<endl;
        }
};

struct curve
{

        VectorXd X;
        VectorXd Y;
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
// extern Generator g;
// extern int Kmax;
// extern int Curve_num;
// extern int Nm;
// extern double w_eps, v0_eps, sigmav2_eps;
// extern ofstream testfile;

/********************************************************************************************************************/
/*Routines. */
/********************************************************************************************************************/
void xixj(MatrixXd & K, const VectorXd &X);
double log_likelihood_micro2(const curve &Dat, const pq_point &theta, int k);
double log_likelihood2(const curve Data[], const pq_point &theta);
double logSumExp(const VectorXd& x);
void draw_initial_model(const curve Data[], pq_point &theta, double *logl);
double prop_split(const curve Data[], pq_point &theta, int k, int *k1, int *k2);
double prop_merge(const curve Data[], pq_point &theta_old, pq_point &theta_new, int *k, int k1, int k2);
double compute_log_split_ratio(pq_point &theta, pq_point &theta_split, int k, int k_split1, int k_split2);
double calc_secondary_moment(const curve &Dat, const pq_point &theta, const pq_point &theta_new, int k, int k1, int k2);
double compute_log_prior_ratio(const pq_point &theta, const pq_point &theta_old);
void zeros(int c[], int K);
