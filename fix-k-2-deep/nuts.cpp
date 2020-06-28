#include "nuts.h"
double log_total_energy(const curve Data[], pq_point &phi, pq_point &theta, const VectorXi &z, double lambda)
{
        int KK = theta.w.size();
        double lte = 0.0;
        for (int k = 0; k < KK; k++)
        {
                lte +=-U(Data, theta, k, z) - K(phi, k, lambda);
        }
        return lte;
}

// Performs one leapfrom step (NUTS paper, Algorithm 1)
void leapfrog(const curve Data[], pq_point &phi, pq_point &theta, double epsilon,
              const VectorXi &z, const double lambda)
{
        //#pragma omp parallel
        //    {
        // theta.print();
        phi.w = phi.w - epsilon*  0.5 * U_partial_w_cpp(Data, theta, z);
        // cout <<"epsilon:"<<epsilon<<endl;

        // phi.print();
        // cout<<"U_partial_w_cpp(Data, theta, z):"<<U_partial_w_cpp(Data, theta, z)<<endl;
        phi.v = phi.v - epsilon * 0.5 * U_partial_v_cpp(Data, theta, z);
        // cout<<"U_partial_v_cpp(Data, theta, z):"<<U_partial_v_cpp(Data, theta, z)<<endl;
        phi.sigma2 = phi.sigma2 - epsilon * 0.5 * U_partial_sigma2_cpp(Data, theta, z);
        // cout<<"U_partial_sigma2_cpp(Data, theta, z):"<<U_partial_sigma2_cpp(Data, theta, z)<<endl;
        theta.w = theta.w + epsilon * phi.w / lambda;
        theta.v = theta.v + epsilon * phi.v / lambda;
        theta.sigma2 = theta.sigma2 + epsilon * phi.sigma2 / lambda;

        theta.positive_reflect();

        phi.w = phi.w - epsilon * 0.5 * U_partial_w_cpp(Data, theta, z);
        phi.v = phi.v - epsilon * 0.5 * U_partial_v_cpp(Data, theta, z);
        phi.sigma2 = phi.sigma2 - epsilon * 0.5 * U_partial_sigma2_cpp(Data, theta, z);
        //}
}

// U-Turn criterion in the generalized form applicable to Riemanian spaces
// See Betancourt's Conceptual HMC (page 58)
bool compute_criterion(const pq_point &p_sharp_minus,
                       const pq_point &p_sharp_plus,
                       const pq_point &rho)
{
        testfile << "p_sharp_plus:" << endl;
        testfile<< p_sharp_plus.w.transpose() << endl;
        testfile << p_sharp_plus.v.transpose() <<endl;
        testfile << p_sharp_plus.sigma2.transpose()<<endl;

        testfile << "p_sharp_minus:" << endl;
        testfile<< p_sharp_minus.w.transpose() << endl;
        testfile << p_sharp_minus.v.transpose() <<endl;
        testfile << p_sharp_minus.sigma2.transpose()<<endl;

        testfile << "rho:" << endl;
        testfile<< rho.w.transpose() << endl;
        testfile <<  rho.v.transpose() <<endl;
        testfile <<  rho.sigma2.transpose()<<endl;

        // p_sharp_plus.write_file("p_sharp_plus", testfile);
        // // testfile<<endl;
        // rho.write_file("rho",testfile);
        // p_sharp_minus.write_file("p_sharp_minus",testfile);
        return (p_sharp_plus % rho) > 0 && (p_sharp_minus % rho) > 0;
}

//**
// Recursively build a new subtree to completion or until the subtree becomes invalid.
// Returns validity of the resulting subtree
// @param z last visited state?
// @param depth Depth of the desired subtree
// @z_propose State proposed from subtree
// @p_sharp_minus p_sharp from left boundary of returned tree (p_sharp = inv(M)*p)
// @p_sharp_plus p_sharp from right boundary of returned tree
// @rho Summed momentum accross trajectory (to compute the generalized stoppin criteria)
int BuildTree(const curve Data[], pq_point &phi, pq_point &theta,
              pq_point &phi_propose, pq_point &theta_propose,
              pq_point &p_sharp_minus,
              pq_point &p_sharp_plus,
              pq_point &rho,
              nuts_util &util,
              int depth, double epsilon,
              const double lambda,
              const VectorXi &z)
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
                util.sum_prob +=  (exp(joint - util.H0)<1) ? exp(joint - util.H0) : 1;
                util.n_tree += 1;
                theta_propose = theta;
                phi_propose = phi;
                rho = rho + phi_propose;
                p_sharp_minus = phi_propose; // p_sharp = inv(M)*p (Betancourt 58)
                p_sharp_plus = p_sharp_minus;
                //  rho.print("rho:");
                return valid_subtree;
        }

        // General recursion
        pq_point p_sharp_dummy;

        // Build the left subtree
        pq_point rho_left(KK);
        rho_left.zeros();
        testfile<<"left tree depth:"<<depth-1<<endl;
        testfile << "p_sharp_minus:" << endl;
        testfile<< p_sharp_minus.w.transpose() << endl;
        testfile << p_sharp_minus.v.transpose() <<endl;
        testfile << p_sharp_minus.sigma2.transpose()<<endl;
        p_sharp_dummy = p_sharp_plus;
        int n1 = BuildTree(Data, phi, theta, phi_propose, theta_propose, p_sharp_minus, p_sharp_plus, rho_left, util, depth - 1, epsilon, lambda, z);
        p_sharp_plus = p_sharp_dummy;
        if (!util.criterion) {
                testfile<<"build exit"<<endl;
                testfile << "p_sharp_plus:" << endl;
                testfile<< p_sharp_plus.w.transpose() << endl;
                testfile << p_sharp_plus.v.transpose() <<endl;
                testfile << p_sharp_plus.sigma2.transpose()<<endl;
                return 0; // early stopping
        }


        // Build the right subtree
        pq_point phi_propose_new(phi_propose);
        pq_point theta_propose_new(theta_propose);
        pq_point rho_right(KK);
        rho_right.zeros();
        testfile<<"right tree depth:"<<depth-1<<endl;
        testfile << "p_sharp_plus:" << endl;
        testfile<< p_sharp_plus.w.transpose() << endl;
        testfile << p_sharp_plus.v.transpose() <<endl;
        testfile << p_sharp_plus.sigma2.transpose()<<endl;
        p_sharp_dummy = p_sharp_minus;
        int n2 = BuildTree(Data, phi, theta, phi_propose_new, theta_propose_new,
                           p_sharp_minus, p_sharp_plus, rho_right, util, depth - 1, epsilon, lambda, z);
        p_sharp_minus = p_sharp_dummy;
        // if (!util.criterion)
        //         return 0;                    // early stopping
        // Choose which subtree to propagate a sample up from.
        //double accept_prob = static_cast<double>(n2) / static_cast<double>(n1 + n2);
        double accept_prob = static_cast<double>(n2) / max((n1 + n2), 1); // avoids 0/0;
        float rand01 = kiss(g);
        if (util.criterion && (rand01 < accept_prob))
        {
                phi_propose = phi_propose_new;
                theta_propose = theta_propose_new;
        }

        // Break when NUTS criterion is no longer satisfied
        //rho_left.print("rho_left:");
        //rho_right.print("rho_right:");
        pq_point rho_subtree = rho_left + rho_right;
        rho = rho + rho_subtree;
        util.criterion = compute_criterion(p_sharp_minus, p_sharp_plus, rho);

        int n_valid_subtree = n1 + n2;
        return (n_valid_subtree);
}
