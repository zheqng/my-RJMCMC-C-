#include "HMC.h"
extern Generator g;
extern int Kmax;
extern int Curve_num,THREAD_NUM;
extern int Nm;
extern double w_eps,v0_eps,sigmav2_eps;
extern ofstream testfile;


struct nuts_util {
        // Constants through each recursion
        double log_u; // uniform sample
        double H0; // Hamiltonian of starting point?
        int sign; // direction of the tree in a given iteration/recursion

        // Aggregators through each recursion
        int n_tree;
        double sum_prob;
        bool criterion;

        // just to guarantee bool initializes to valid value
        nuts_util() : criterion(false) {
        }
};

double log_total_energy(const curve Data[],pq_point &phi, pq_point &theta,const VectorXi & z,double lambda);
void leapfrog(const curve Data[], pq_point &phi, pq_point &theta, double epsilon,
              const VectorXi & z, const double lambda);
bool compute_criterion(const pq_point & p_sharp_minus,
                       const pq_point & p_sharp_plus,
                       const pq_point & rho);
int BuildTree(const curve Data[], pq_point& phi, pq_point &theta,
              pq_point& phi_propose,pq_point &theta_propose,
              pq_point & p_sharp_left,
              pq_point & p_sharp_right,
              pq_point & rho,
              nuts_util& util,
              int depth, double epsilon,
              const double lambda,
              const VectorXi & z);
void  sample_nuts_cpp(const curve Data[], pq_point & current_theta,
                      const VectorXi & z);
// void find_reasonable_epsilon(const curve Data[],pq_point & phi,
//   pq_point &theta,double & epsilon,const vec & z,const double lambda);
