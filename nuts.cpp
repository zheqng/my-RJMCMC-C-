#include "nuts.h"
double log_total_energy(const curve Data[], pq_point &phi, pq_point &theta, const vec &z, double lambda)
{
        int KK = theta.w.size();
        double lte = 0.0;
        for (int k = 0; k < KK; k++)
        {
                lte +=-U(Data, theta, k, z) + K(phi, k, lambda);
        }
        return lte;
}

// Performs one leapfrom step (NUTS paper, Algorithm 1)
void leapfrog(const curve Data[], pq_point &phi, pq_point &theta, double epsilon,
              const vec &z, const double lambda)
{
        phi.w = phi.w - epsilon * 0.5 * U_partial_w_cpp(Data, theta, z);
        phi.v = phi.v - epsilon * 0.5 * U_partial_v_cpp(Data, theta, z);
        phi.sigma2 = phi.sigma2 - epsilon * 0.5 * U_partial_sigma2_cpp(Data, theta, z);

        theta.w = theta.w + epsilon * phi.w / lambda;
        theta.v = theta.v + epsilon * phi.v / lambda;
        theta.sigma2 = theta.sigma2 + epsilon * phi.sigma2 / lambda;

        theta.positive_reflect();

        phi.w = phi.w - epsilon * 0.5 * U_partial_w_cpp(Data, theta, z);
        phi.v = phi.v - epsilon * 0.5 * U_partial_v_cpp(Data, theta, z);
        phi.sigma2 = phi.sigma2 - epsilon * 0.5 * U_partial_sigma2_cpp(Data, theta, z);

}

// U-Turn criterion in the generalized form applicable to Riemanian spaces
// See Betancourt's Conceptual HMC (page 58)
bool compute_criterion(const pq_point &p_sharp_minus,
                       const pq_point &p_sharp_plus,
                       const pq_point &rho)
{
        return (p_sharp_plus % rho).accu() > 0 && (p_sharp_minus % rho).accu() > 0;
}

//**
// Recursively build a new subtree to completion or until the subtree becomes invalid.
// Returns validity of the resulting subtree
// @param z last visited state?
// @param depth Depth of the desired subtree
// @z_propose State proposed from subtree
// @p_sharp_left p_sharp from left boundary of returned tree (p_sharp = inv(M)*p)
// @p_sharp_right p_sharp from right boundary of returned tree
// @rho Summed momentum accross trajectory (to compute the generalized stoppin criteria)
int BuildTree(const curve Data[], pq_point &phi, pq_point &theta,
              pq_point &phi_propose, pq_point &theta_propose,
              pq_point &p_sharp_left,
              pq_point &p_sharp_right,
              pq_point &rho,
              nuts_util &util,
              int depth, double epsilon,
              const double lambda,
              const vec &z)
{

        //cout << "\n Tree direction:" << util.sign << " Depth:" << depth << endl;

        int KK = theta.w.size();


        float delta_max = 1000; // Recommended in the NUTS paper: 1000

        // Base case - take a single leapfrog step in the direction v
        if (depth == 0)
        {

                leapfrog(Data, phi, theta, (util.sign) * epsilon, z, lambda);

                float joint = log_total_energy(Data, phi, theta, z, lambda);
                int valid_subtree = (util.log_u <= joint); // Is the new point in the slice?
                util.criterion = util.log_u - joint < delta_max; // Is the simulation wildly inaccurate? // TODO: review
                util.sum_prob +=  (exp(-joint + util.H0)<1) ? exp(joint - util.H0) : 1;
                util.n_tree += 1;
                theta_propose = theta;
                phi_propose = phi;
                rho = rho + phi;
                p_sharp_left = phi; // p_sharp = inv(M)*p (Betancourt 58)
                p_sharp_right = p_sharp_left;
                return valid_subtree;
        }

        // General recursion
        pq_point p_sharp_dummy;

        // Build the left subtree
        pq_point rho_left(KK);
        rho_left.zeros();
        int n1 = BuildTree(Data, phi, theta, phi_propose, theta_propose, p_sharp_left, p_sharp_dummy, rho_left, util, depth - 1, epsilon, lambda, z);
        if (!util.criterion)
                return 0; // early stopping

        // Build the right subtree
        pq_point phi_propose_right(phi);
        pq_point theta_propose_right(theta);
        pq_point rho_right(KK);
        rho_left.zeros();
        int n2 = BuildTree(Data, phi, theta, phi_propose_right, theta_propose_right,
                           p_sharp_dummy, p_sharp_right, rho_right, util, depth - 1, epsilon, lambda, z);
        // Choose which subtree to propagate a sample up from.
        //double accept_prob = static_cast<double>(n2) / static_cast<double>(n1 + n2);
        double accept_prob = static_cast<double>(n2) / max((n1 + n2), 1); // avoids 0/0;
        float rand01 = kiss(g);
        if (util.criterion && (rand01 < accept_prob))
        {
                phi_propose = phi_propose_right;
                theta_propose = theta_propose_right;
        }

        // Break when NUTS criterion is no longer satisfied
        pq_point rho_subtree = rho_left + rho_right;
        rho = rho + rho_subtree;
        util.criterion = compute_criterion(p_sharp_left, p_sharp_right, rho);

        int n_valid_subtree = n1 + n2;
        return (n_valid_subtree);
}

void sample_nuts_cpp(const curve Data[], pq_point &current_theta,
                     const vec &z)
{


        int iter = 2;
        int MAXDEPTH = 5;
        double lambda = 1.0;
        int KK = current_theta.w.size();

        int I_adapt = 10, t0 = iter;
        double delta = 0.65, kappa = 0.75, gamma = 0.05;
        double H_bar[100], log_eps_bar[100], log_eps[100];
        double epsilon = 1e-6;
        log_eps[0] = log(epsilon);
        H_bar[0] = 0;
        double mu = log(10.0 * epsilon);


        // Store fixed data and parameters

        pq_point samples[iter];
        pq_point phi_0(KK); // initial momentum
        // current_theta = log(current_theta); // Transform to unrestricted space
        samples[0] = current_theta;

        pq_point phi;
        pq_point theta;

        // Used to compute the NUTS generalized stopping criterion
        pq_point rho(KK);
        nuts_util util;

        // Transition
        for (int i = 1; i < iter; i++)
        {
                cout << "******************* Sample: " << i << endl;

                // Sample new momentum (K independent standard normal variates)
                phi_0.initialize(g);

                // Initialize the path. Proposed sample,
                // and leftmost/rightmost position and momentum
                ////////////////////////
                theta = current_theta;
                phi = phi_0;
                pq_point phi_plus(phi);
                pq_point theta_plus(theta);
                pq_point phi_minus(phi);
                pq_point theta_minus(theta);
                pq_point phi_propose(phi);
                pq_point theta_propose(theta);

                // Utils o compute NUTS stop criterion
                pq_point p_sharp_plus = phi;
                pq_point p_sharp_dummy = p_sharp_plus;
                pq_point p_sharp_minus = p_sharp_plus;
                pq_point rho(phi);

                // Hamiltonian
                // Joint logprobability of position q and momentum p
                float joint = log_total_energy(Data, phi_0, current_theta, z, lambda);
                util.H0 = joint;

                // Slice variable
                ///////////////////////
                // Sample the slice variable: u ~ uniform([0, exp(joint)]).
                // Equivalent to: (log(u) - joint) ~ exponential(1).
                // logu = joint - exprnd(1);
                float random = -log(kiss(g));
                util.log_u = joint - random;

                int n_valid = 1;
                util.criterion = true;

                // Build a trajectory until the NUTS criterion is no longer satisfied
                int depth_ = 0;
                int divergent_ = 0;
                util.n_tree = 0;
                util.sum_prob = 0;

                // Build a balanced binary tree until the NUTS criterion fails
                while (util.criterion && (depth_ < MAXDEPTH))
                {

                        // Build a new subtree in the chosen direction
                        // (Modifies z_propose, z_minus, z_plus)
                        pq_point rho_subtree(KK);
                        rho_subtree.zeros();

                        // Build a new subtree in a random direction
                        util.sign = 2 * (kiss(g) < 0.5) - 1;
                        int n_valid_subtree = 0;
                        if (util.sign == 1)
                        {
                                phi.pq_point::operator=(phi_minus);
                                theta.pq_point::operator=(theta_minus);
                                n_valid_subtree = BuildTree(Data, phi, theta, phi_propose, theta_propose,
                                                            p_sharp_dummy, p_sharp_plus, rho_subtree, util, depth_, epsilon, lambda, z);
                                phi_plus.pq_point::operator=(phi);
                                theta_plus.pq_point::operator=(theta);
                        }
                        else
                        {
                                phi = (phi_plus);
                                theta = (theta_plus);
                                n_valid_subtree = BuildTree(Data, phi, theta, phi_propose, theta_propose,
                                                            p_sharp_dummy, p_sharp_minus, rho_subtree, util, depth_, epsilon, lambda, z);
                                phi_minus.pq_point::operator=(phi);
                                theta_minus.pq_point::operator=(theta);
                        }
                        //if(!valid_subtree) break;
                        ++depth_; // Increment depth.
                        if (util.criterion)
                        {
                                // Use Metropolis-Hastings to decide whether or not to move to a
                                // point from the half-tree we just generated.
                                double subtree_prob = min(1.0, static_cast<double>(n_valid_subtree) / n_valid);
                                if (kiss(g) < subtree_prob)
                                {
                                        current_theta = theta_propose; // Accept proposal (it will be THE new sample when s=0)
                                }
                        }

                        // Update number of valid points we've seen.
                        n_valid += n_valid_subtree;

                        // Break when NUTS criterion is no longer satisfied
                        rho = rho + rho_subtree;
                        util.criterion = util.criterion && compute_criterion(p_sharp_minus, p_sharp_plus, rho);
                } // end while

                samples[i] = current_theta;
                if (i < I_adapt)
                {
                        H_bar[i] = (1.0 - 1.0 / (i + t0)) * H_bar[i - 1] + 1.0 / (i + t0) * (delta - (util.sum_prob) / (util.n_tree));
                        // testfile<<H_bar[i]<<" ";
                        log_eps[i] = mu - sqrt((double)i) / gamma * H_bar[i];
                        epsilon = exp(log_eps[i]);
                        log_eps_bar[i] = pow(i, -kappa) * log_eps[i] + (1 - pow(i, -kappa)) * log_eps_bar[i - 1];
                }
                else
                {
                        log_eps[i] = log_eps_bar[I_adapt];
                        epsilon = exp(log_eps[i]);
                }
                // testfile << (util.sum_prob) / (util.n_tree) << " " << epsilon << endl;
        } // end for

        // return samples;
}
