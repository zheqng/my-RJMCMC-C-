/********************************************************************************************************************/
/*mix_lib: Library of functions used by rj_mix. */
/********************************************************************************************************************/
#include "mix_lib.h"
void xixj(MatrixXd &K, const VectorXd &X)
{
								int r = X.size();
								// mat K(r, r);
								// #pragma omp parallel for
								for (int i = 0; i < r; ++i)
								{
																for (int j = i; j < r; ++j)
																{
																								K(j, i) = pow(X(j) - X(i), 2.0);
																								K(i, j) = K(j, i);
																}
								}
}

double log_likelihood_micro2(const curve &Dat, const pq_point &theta, int k)
{


								int r = Dat.X.size();
								MatrixXd cov(r, r), SquareX(r, r);

// caculate cov
								MatrixXd EYE = MatrixXd::Identity(r, r);
								xixj(SquareX, Dat.X);
								MatrixXd Inner = -SquareX * theta.w(k) / 2.0;
								cov = theta.v(k) * Inner.array().exp().matrix() + EYE * theta.sigma2(k);
// testfile << cov << endl;
								double tmp;
// testfile<<(cov * solve(cov, Dat.Y))<<endl;
								LLT<MatrixXd> llt;
								llt.compute(cov);
// double tmp =;
								tmp = -Dat.Y.dot(llt.solve(Dat.Y))   / 2.0 - (double)r / 2.0 * log(2.0 * PI) - 1.0 / 2.0 * log(abs(cov.determinant()));
// cout<<"tmp1:"<<Dat.Y.dot(llt.solve(Dat.Y))<<" tmp3:"<<log(abs(cov.determinant()))<<endl;
// cout<<"tmp:"<<tmp<<endl;
								return tmp;
}
// more than one operator "/" matches these operands:
//  -- function template "arma::enable_if2<arma::is_arma_type<T1>::value,
//  const arma::eOp<T1, arma::eop_scalar_div_post>>::result arma::operator/(const T1 &X, T1::elem_type k)"
// -- function template "arma::enable_if2<arma::is_arma_type<T1>::value,
// const arma::eOp<T1, arma::eop_scalar_div_post>>::result operator/(const T1 &X, T1::elem_type k)"
// -- operand types are: const arma::eOp<arma::eOp<arma::mat, arma::eop_neg>, arma::eop_scalar_times> / double

// double log_likelihood2(const curve Data[], const pq_point &theta)
// {
//         int t, k;
//         int K = theta.w.size();
//         double norm;
//         // vec logP(K);
//         norm = 0.0;
// // time_t t_start, t_end;
// // // ofstream TimeFP("time.txt");
// // double DiffTime;
// // t_start = time(NULL);
// // pq_poinSt theta;
//         omp_set_nested(1);
//         // #pragma omp parallel for reduction (+:norm)
//         for (t = 0; t < Curve_num; t++)
//         {
//                 // cout<<"The outer loop:thread_num"<<omp_get_thread_num()<<endl;
//                 vec logP = zeros<vec>(K);
//                 // #pragma omp parallel for
//                 for (k = 0; k < K; k++)
//                 {
//                         // logP(k) = log_likelihood_micro2(Data[t], theta, k) + log(theta.pi(k));
//                         logP(k) = log_likelihood_micro2(Data[t], theta, k);
//
//                 }
//                 logP = logP +  log(theta.pi);
//                 norm = norm + logSumExp(logP);
//                 // printf("ID: %d, Max threads: %d, Num threads: %d \n",omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());
//
//         }
//         // t_end = time(NULL);
//         // DiffTime = difftime(t_end, t_start);
//         // cout<<DiffTime<<endl;
//         return norm;
// }

double log_likelihood2(const curve Data[], const pq_point &theta)
{
								// int t, k;
								int K = theta.w.size();
								double norm;
								// vec logP(K);
								norm = 0.0;
// time_t t_start, t_end;
// // ofstream TimeFP("time.txt");
// double DiffTime;
// t_start = time(NULL);
// pq_poinSt theta;
								omp_set_nested(1);
								#pragma omp parallel for default(none)  reduction(+:norm) shared(Data,theta,K,Curve_num)
								for (int t = 0; t < Curve_num; t++)
								{
																// cout<<"The outer loop:thread_num"<<omp_get_thread_num()<<endl;
																VectorXd logP = VectorXd::Zero(K);

																#pragma omp parallel for default(none) shared(logP,Data,theta,t,K) ordered
																for (int k = 0; k < K; k++)
																{
																								// cout<<"the inner loop: thread_num"<<omp_get_thread_num()<<endl;

																								logP(k) = log_likelihood_micro2(Data[t], theta, k);
																}
																logP = logP + theta.pi.array().log().matrix();;
																norm += logSumExp(logP);
																// printf("ID: %d, Max threads: %d, Num threads: %d \n",omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());

								}
								// t_end = time(NULL);
								// DiffTime = difftime(t_end, t_start);
								// cout<<DiffTime<<endl;
								return norm;
}

double logSumExp(const VectorXd &x)
{
								LDOUBLE maxv = x.maxCoeff();

								LDOUBLE cumsum = 0.0;
								for (int i = 0; i < x.size(); i++)
																cumsum += exp(x(i) - maxv);

								return maxv + log(cumsum);
}

/************************************************************************/
/* Initializes the parameters of the model (and returns the scaling	*/
/* needed  to obtain a non -Inf log-likelihood, this should normaly be	*/
/* 1)									*/
/************************************************************************/

void draw_initial_model(const curve Data[], pq_point &theta, double *logl)
{

								/* Start by drawing from prior					*/

								theta.draw_from_prior(g);

								*logl = log_likelihood2(Data, theta);
}

/************************************************************************/
/* Proposes a split move and returns log-likelihood			*/
/* put new classes in current and last position*/
/************************************************************************/
double prop_split(const curve Data[], pq_point &theta, int k, int *k1, int *k2)
{
								double u[4];
								double logl;
								pq_point Elem_mix1(1), Elem_mix2(1);
								int K = theta.w.size();

								if (K >= Curve_num)
								{

																exit(1);
								}
								if ((k < 0) || (k >= K))
								{
																exit(1);
								}

								u[0] = Beta(2.0, 2.0, g);
								u[1] = Beta(1.0, 1.0, g);
								u[2] = Beta(2.0, 2.0, g);
								u[3] = Beta(2.0, 2.0, g);

								/*_____________________________reallocate theta___________________*/
								double pi_k = theta.pi(k);
								Elem_mix1.pi(0) = u[0] * theta.pi(k);
								Elem_mix2.pi(0) = (1.0 - u[0]) * theta.pi(k);

								Elem_mix1.sigma2(0) = u[1] * theta.sigma2(k) * pi_k / Elem_mix1.pi(0);
								Elem_mix2.sigma2(0) = (1.0 - u[1]) * theta.sigma2(k) * pi_k / Elem_mix2.pi(0);

								Elem_mix1.v(0) = u[2] * theta.v(k) * pi_k / Elem_mix1.pi(0);
								Elem_mix2.v(0) = (1.0 - u[2]) * theta.v(k) * pi_k / Elem_mix2.pi(0);

								Elem_mix1.w(0) = (1.0 - u[3]) / u[2] * theta.w(k);
								Elem_mix2.w(0) = u[3] / (1.0 - u[2]) * theta.w(k);

								/*____________________________insert w[k]________________________*/
								// m.deleteP_seq(k);
								*k1 = theta.w.size();
								*k2 = (*k1) + 1;
								// cout<<"theta old:"<<theta.v<<endl;
								// cout<<"insert:"<<Elem_mix1.v(0)<<" "<<Elem_mix2.v(0)<<endl;
								// cout<<"k:"<<k<<" k1:"<<*k1<<" k2:"<<*k2<<endl;
								theta.insertPre_seq(*k1, Elem_mix1);
								// cout<<"theta new:"<<theta.v<<endl;
								theta.insertPre_seq(*k2, Elem_mix2);
								// cout<<"theta new:"<<theta.v<<endl;
								theta.deleteP_seq(k);
								// cout<<"after delete k-th,theta new:"<<theta.v<<endl;
								*k1 = (*k1) - 1;
								*k2 = (*k2) - 1;
								// cout<<"k:"<<k<<" k1:"<<*k1<<" k2:"<<*k2<<endl;
								theta.pi = theta.pi/theta.pi.sum();

								logl = log_likelihood2(Data, theta);
								return logl;
}

/************************************************************************/
/* Proposes a merge move and returns log-likelihood for new parameters  */
/* (puts the merged class in last position)				*/
/************************************************************************/
double prop_merge(const curve Data[], pq_point &theta_old, pq_point &theta_new, int *k, int k1, int k2)
{
								int i_old, i_new;
								double logl;
								pq_point Elem_mix(1);
								if (theta_old.w.size() <= 1)
								{
																cout<<"theta_old.w.size() <= 1;"<<endl;
																exit(1);
								}
								if ((k1 < 0) || (k2 < 0) || (k1 >= theta_old.w.size()) || (k2 >= theta_old.w.size()) || (k1 == k2))
								{
																cout<<"(k1 < 0) || (k2 < 0) || (k1 >= theta_old.w.size()) || (k2 >= theta_old.w.size()) || (k1 == k2)"<<endl;
																exit(1);
								}

								/* Merge classes k1 and k2 in the last one of the new model		*/
								Elem_mix.pi(0) = theta_old.pi(k1) + theta_old.pi(k2);

								Elem_mix.sigma2(0) = (theta_old.pi(k1) * theta_old.sigma2(k1) +
																														theta_old.pi(k2) * theta_old.sigma2(k2)) /
																													Elem_mix.pi(0);

								Elem_mix.v(0) = (theta_old.pi(k1) * theta_old.v(k1) +
																									theta_old.pi(k2) * theta_old.v(k2)) /
																								Elem_mix.pi(0);
								// cout<<"theta old:"<<theta_old.v<<endl;
								// cout<<"theta new:"<<Elem_mix.v(0)<<endl;
								Elem_mix.w(0) = (theta_old.pi(k1) * theta_old.v(k1) * theta_old.w(k1) +
																									theta_old.pi(k2) * theta_old.v(k2) * theta_old.w(k2)) /
																								Elem_mix.pi(0) / Elem_mix.v(0);

								/*_______________________merge_________________________*/
								theta_new = theta_old;

								*k = theta_new.w.size();

								theta_new.insertPre_seq(*k, Elem_mix);
								// cout<<"theta new:"<<theta_new.v<<endl;

								//delete k1 and k2
								theta_new.deleteP_seq(k1);
								theta_new.deleteP_seq(k2-1);
								// cout<<"theta new:"<<theta_new.v<<endl;
								*k = (*k) - 2;
								theta_new.pi = theta_new.pi/theta_new.pi.sum();
								// cout<<"K:"<<*k<<endl;
								// theta_new.print("theta_new");
								/* Compute new log-likelihood						*/
								logl = log_likelihood2(Data, theta_new);
								// cout<<1<<endl;
								return logl;
}

/************************************************************************/
/* Computes the part of the split ratio */
/************************************************************************/
double compute_log_split_ratio(pq_point &theta, pq_point &theta_split, int k, int k_split1, int k_split2)
{
								double lograt = 0.0;
								double u[4];
								int K = theta.w.size();
								// cout<<5<<endl;
								/*_________________________ Weights,pi__________________________*/
								u[0] = theta_split.pi(k_split1) / theta.pi(k);

								lograt += log((double)K) - log(6.0) - log(u[0] * (1 - u[0])) + log(theta.pi(k));
								// cout<<6<<endl;
								/*________________________ sigmav2s_____________________________*/
								u[1] = theta_split.sigma2(k_split1) * theta_split.pi(k_split1) / theta.sigma2(k) / theta.pi(k);
								boost::math::inverse_gamma_distribution<> IG_dist_s(10, 0.01);
								lograt += log(pdf(IG_dist_s, theta_split.sigma2(k_split1))) + log(pdf(IG_dist_s, theta_split.sigma2(k_split2))) - log(pdf(IG_dist_s, theta.sigma2(k)))

								          //-log(6.0) - log(u[1] * (1 - u[1]))
																		+ log(theta.sigma2(k)) + 2.0 * log(theta.pi(k)) - log(theta_split.pi(k_split1) * theta_split.pi(k_split2));
								// cout<<7<<endl;
								/*________________________v0s________________________________*/
								u[2] = theta_split.v(k_split1) * theta_split.pi(k_split1) / theta.v(k) / theta.pi(k);
								boost::math::inverse_gamma_distribution<> IG_dist_v(0.5, 0.5);
								lograt += log(pdf(IG_dist_v, theta_split.v(k_split1))) + log(pdf(IG_dist_v, theta_split.v(k_split2))) - log(pdf(IG_dist_v, theta.v(k)))

																		-log(6.0) - log(u[2] * (1 - u[2])) + log(theta.v(k)) + 2.0 * log(theta.pi(k)) - log(theta_split.pi(k_split1) * theta_split.pi(k_split2));
								// cout<<8<<endl;
								/*________________________w________________________________*/
								u[3] = theta_split.w(k_split2) / theta.w(k) * (1 - u[2]);
								boost::math::inverse_gamma_distribution<> IG_dist_w(0.5, 0.5);
								lograt += log(pdf(IG_dist_w, theta_split.w(k_split1))) + log(pdf(IG_dist_w, theta_split.w(k_split2))) - log(pdf(IG_dist_w, theta.w(k)))

																		-log(6.0) - log(u[3] * (1 - u[3])) + log(theta.w(k)) - log(u[2] * (1 - u[2]));

								return lograt;
}

double calc_secondary_moment(const curve &Dat, const pq_point &theta, const pq_point &theta_new, int k, int k1, int k2)
{
								double secondary_moment = 0.0;
// evaluate the covariance for the mth batch/curve
// there are totally Nm points on this batch/curve
								int r = Nm;
								MatrixXd cov(r, r), cov1(r, r), cov2(r, r), DiffCov(r, r), SquareX(r, r);

// vec colones = ones<vec>(r);
// rowvec rowones = ones<rowvec>(r);
								MatrixXd EYE = MatrixXd::Identity(r, r);

								xixj(SquareX, Dat.X);
//  square(Dat.X * rowones - colones * trans(Dat.X));
								MatrixXd Inner = -SquareX * theta.w(k) / 2.0;
								cov = theta.v(k) * Inner.array().exp().matrix();
								cov = cov + EYE * theta.sigma2(k);
								MatrixXd Inner1 = -SquareX * theta_new.w(k1) / 2.0;
								cov1 = theta_new.v(k1) * Inner.array().exp().matrix() + EYE * theta_new.sigma2(k1);
								MatrixXd Inner2 = -SquareX * theta_new.w(k2) / 2.0;
								cov2 = theta_new.v(k2) * Inner.array().exp().matrix() + EYE * theta_new.sigma2(k2);

								DiffCov = theta.pi(k) * cov - theta_new.pi(k1) * cov1 - theta_new.pi(k2) * cov2;

								return secondary_moment = DiffCov.norm();
}

double compute_log_prior_ratio(const pq_point &theta, const pq_point &theta_old)
{
								double lograt = 0.0;
								int old_K = theta_old.w.size();
								int K = theta.w.size();

								/*_________________________ Weights,pi__________________________*/
								if (old_K < K)
																lograt -= log((double)old_K);
								if (old_K > K)
																lograt += log((double)K);
								if (old_K == K)
																lograt += 0.0;

								/*________________________ sigmav2s_____________________________*/
								boost::math::inverse_gamma_distribution<> IG_dist_s(10, 0.01);
								for (int i = 0; i < old_K; i++)
								{
																lograt = lograt - log(pdf(IG_dist_s, theta_old.sigma2(i)));
								}
								for (int i = 0; i < K; i++)
								{
																lograt = lograt + log(pdf(IG_dist_s, theta.sigma2(i)));
								}

								/*________________________v0s________________________________*/
								boost::math::inverse_gamma_distribution<> IG_dist_v(0.5, 0.5);
								for (int i = 0; i < old_K; i++)
								{
																lograt = lograt - log(pdf(IG_dist_v, theta_old.v(i)));
								}
								for (int i = 0; i < K; i++)
								{
																lograt = lograt + log(pdf(IG_dist_v, theta.v(i)));
								}

								/*________________________w________________________________*/
								boost::math::inverse_gamma_distribution<> IG_dist(0.5, 0.5);
								for (int i = 0; i < old_K; i++)
								{
																lograt = lograt - log(pdf(IG_dist, theta_old.w(i)));
								}
								for (int i = 0; i < K; i++)
								{
																lograt = lograt + log(pdf(IG_dist, theta.w(i)));
								}

								return lograt;
}


void zeros(int c[], int K)
{
								for (int k = 0; k < K; k++)
								{
																c[k] = 0;
								}
}
