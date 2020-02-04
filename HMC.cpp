
#include "HMC.h"

double K(const pq_point &theta, int k, const double lambda)
{
        double kinetic_energy = 0.0;
        kinetic_energy += pow(theta.w(k), 2.0) / 2.0 / lambda + pow(theta.v(k), 2.0) / 2.0 / lambda + pow(theta.sigma2(k), 2.0) / 2.0 / lambda;
        return kinetic_energy;
}

double U(const curve Data[], const pq_point &theta, int k, const VectorXi &z)
{
        double u = 0.0, l_pk = 0.0;
        int zk[MaxM];
        int t=0;
        for(int m=0; m<Curve_num; m++) {
                if(z(m)==k) {zk[t]=m; t++;}
        }

        // VectorXd l_pk = VectorXd::Zero(t);
        int m;
        #pragma omp parallel for  default(none) ordered reduction(+:l_pk) shared(t,Data,zk,theta,k) private(m)
        for ( m = 0; m < t; m++)
        {
                l_pk +=   log_likelihood_micro2(Data[zk[m]], theta, k);
        }


        boost::math::inverse_gamma_distribution<> IG_dist_w(0.5, 0.5);
        boost::math::inverse_gamma_distribution<> IG_dist_v(0.5, 0.5);
        boost::math::inverse_gamma_distribution<> IG_dist_s(10, 0.01);

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
        MatrixXd cov(r, r), covinv(r, r), cov_partial_w(r, r), squareX(r, r);
        double tmp;
        VectorXd covinv_Y(r);

        // vec colones = ones<vec>(r);
        // rowvec rowones = ones<rowvec>(r);
        MatrixXd EYE = MatrixXd::Identity(r, r);
        xixj(squareX,Dat.X);
        // square(Dat.X * rowones - colones * trans(Dat.X));
        cov = theta.v(k) * ((-squareX * theta.w(k) / 2.0).array().exp().matrix()) + EYE * theta.sigma2(k);
        // covinv = inv(cov);
        cov_partial_w = -0.5 * theta.v(k) * ( (squareX.array()) * ((-squareX * theta.w(k) / 2.0).array().exp()) ).matrix();
        // tmp = -0.5 * trace(covinv * cov_partial_w);
        LLT<MatrixXd> llt;
        llt.compute(cov);
        tmp=  -0.5*llt.solve(cov_partial_w).trace();
        // tmp = -0.5 * trace(solve(cov, cov_partial_w));
        covinv_Y = llt.solve(Dat.Y);
        // tmp = tmp + 0.5 * trans(Dat.Y) * covinv * cov_partial_w * covinv * Dat.Y;
        tmp = tmp + 0.5 * covinv_Y.transpose() * cov_partial_w * covinv_Y;
        return tmp;
}

VectorXd U_partial_w_cpp(const curve Data[], const pq_point &theta,  const VectorXi &z)
{
        int K = theta.w.size();
        VectorXd u_partial_w =  VectorXd::Zero(K);
        int zk[MaxM];
        int t=0;

        omp_set_nested(1);

        #pragma omp parallel for default(none) shared(z,u_partial_w,Data,theta,K,Curve_num) private(t,zk)
        for (int k = 0; k < K; k++)
        {
                t=0;
                #pragma omp parallel for  default(none) ordered shared(z,zk,t,k,Curve_num)
                for(int m=0; m<Curve_num; m++) {
                        if(z(m)==k) {zk[t]=m; t++;}
                }

                #pragma omp parallel for  default(none) ordered  shared(u_partial_w,t,Data,zk,theta,k)
                for (int m = 0; m < t; m++)
                {
                        u_partial_w(k) = u_partial_w(k) - logP_partial_w_percurve_cpp(Data[zk[m]], theta, k);
                }
                // cout << " u_partial_w:"<<u_partial_w<<endl;
        }
        u_partial_w += (0.5 + 1) * theta.w.cwiseInverse() - 0.5 * (theta.w.array().square().matrix().cwiseInverse());
        // cout << " u_partial_w:"<<u_partial_w<<endl;
        // u_partial_w + (0.5 + 1) / theta.w - 0.5 / square(theta.w);

        return u_partial_w;
}
/*U0 partial v0*/
double logP_partial_v_percurve_cpp(const curve &Dat, const pq_point &theta, int k)
{
        int r = Dat.X.size();
        MatrixXd cov(r, r), covinv(r, r), cov_partial_v(r, r), squareX(r, r);
        double tmp;
        VectorXd covinv_Y(r);

        // vec colones = ones<vec>(r);
        // rowvec rowones = ones<rowvec>(r);
        MatrixXd EYE = MatrixXd::Identity(r, r);
        xixj(squareX,Dat.X);
        // square(Dat.X * rowones - colones * trans(Dat.X));

        cov = theta.v(k) * ((-squareX * theta.w(k) / 2.0).array().exp().matrix() ) + EYE * theta.sigma2(k);
        // covinv = inv(cov);
        cov_partial_v = (-squareX * theta.w(k) / 2.0).array().exp().matrix();
        LLT<MatrixXd> llt;
        llt.compute(cov);
        covinv_Y = llt.solve(Dat.Y);
        tmp = -0.5 * llt.solve(cov_partial_v).trace() + 0.5 * covinv_Y.transpose() * cov_partial_v * covinv_Y;

        return tmp;
}

VectorXd U_partial_v_cpp(const curve Data[], const pq_point &theta, const VectorXi &z)
{
        int K = theta.w.size();
        VectorXd u_partial_v =  VectorXd::Zero(K);
        omp_set_nested(1);
        int zk[MaxM];
        int t;
        #pragma omp parallel for default(none) shared(z,u_partial_v,Data,theta,K,Curve_num) private(t,zk)
        for (int k = 0; k < K; k++)
        {
                t=0;
                #pragma omp parallel for  default(none) ordered shared(z,zk,t,k,Curve_num)
                for(int m=0; m<Curve_num; m++) {
                        if(z(m)==k) {zk[t]=m; t++;}
                }
                #pragma omp parallel for  default(none) ordered  shared(u_partial_v,t,Data,zk,theta,k)
                for (int m = 0; m < t; m++)
                {
                        u_partial_v(k) = u_partial_v(k) - logP_partial_v_percurve_cpp(Data[zk[m]], theta, k);
                }

        }
        u_partial_v+= (0.5 + 1) * theta.v.cwiseInverse() - 0.5 * (theta.v.array().square().matrix().cwiseInverse());
        return u_partial_v;
}


/*U0 partial sigmav2*/
double logP_partial_sigma2_percurve_cpp(const curve &Dat, const pq_point &theta, int k)
{
        int r = Dat.X.size();
        MatrixXd cov(r, r), covinv(r, r), cov_partial_sigma2(r, r), squareX(r, r);
        double tmp;
        VectorXd covinv_Y(r);

        MatrixXd EYE =  MatrixXd::Identity(r, r);
        xixj(squareX,Dat.X);

        cov = theta.v(k) * ((-squareX * theta.w(k) / 2.0).array().exp().matrix()) + EYE * theta.sigma2(k);
// covinv = inv(cov);
        cov_partial_sigma2 = EYE;
        LLT<MatrixXd> llt;
        llt.compute(cov);
        covinv_Y = llt.solve(Dat.Y);
        tmp = -0.5 * llt.solve(cov_partial_sigma2).trace() + 0.5 * covinv_Y.transpose() *  cov_partial_sigma2 * covinv_Y;

        return tmp;
}

VectorXd U_partial_sigma2_cpp(const curve Data[], const pq_point &theta, const VectorXi &z)
{
        int K = theta.w.size();
        VectorXd u_partial_sigma2 = VectorXd::Zero(K);
        int zk[MaxM];
        int t;
        omp_set_nested(1);

        #pragma omp parallel for default(none) ordered shared(z,u_partial_sigma2,Data,theta,K,Curve_num) private(t,zk)
        for (int k = 0; k < K; k++)
        {
                t=0;
                #pragma omp parallel for  default(none) ordered shared(z,zk,t,k,Curve_num)
                for(int m=0; m<Curve_num; m++) {
                        if(z(m)==k) {zk[t]=m; t++;}
                }

                #pragma omp parallel for  default(none) ordered  shared(u_partial_sigma2,t,Data,zk,theta,k)
                for (int m = 0; m < t; m++)
                {
                        u_partial_sigma2(k) = u_partial_sigma2(k) - logP_partial_sigma2_percurve_cpp(Data[zk[m]], theta, k);
                }
        }



        u_partial_sigma2   += (10.0 + 1.0) * theta.sigma2.cwiseInverse() - 0.01 * (theta.sigma2.array().square().matrix().cwiseInverse());

        return u_partial_sigma2;
}

/*HMC */

/*______________________________Gibbs Sampling____________________________*/

void Gibbs_Sampling_pi(pq_point &theta, const VectorXi &z)
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

        theta.pi = theta.pi/theta.pi.sum();
}
int max_index(VectorXd logp_k){
        double max_value;
        int max_id=0;
        max_value = logp_k(0);
        for(int i=1; i<logp_k.size(); i++) {
                if(max_value<logp_k(i)) {max_id = i; max_value = logp_k(i);}
        }
        return max_id;
}

void Gibbs_Sampling_z(const curve Data[], const pq_point &theta, VectorXi &z)
{
        int K = theta.w.size();
        // vec logp_k(K);
        z = VectorXi::Zero(Curve_num);
        omp_set_nested(1);

        #pragma omp parallel for  default(none) shared(Curve_num,K,Data,theta,z)
        for (int m = 0; m < Curve_num; m++)
        {
                VectorXd logp_k(K);
                #pragma omp parallel for default(none) ordered shared(K,Data,theta,logp_k,m)
                for (int k = 0; k < K; k++)
                {
                        logp_k(k) =log_likelihood_micro2(Data[m], theta, k);
                }
                /*_________________________hard cut___________________________*/
                logp_k= logp_k +  theta.pi.array().log().matrix();
                z(m) = max_index(logp_k);
        }
}
