/********************************************************************************************************************/
/*mix_lib: Library of functions used by rj_mix. */
/********************************************************************************************************************/
#include "mix_lib.h"
void  xixj(mat & K,const vec &X)
{
	int r = X.size();
	// mat K(r, r);
	for (int i = 0; i < r; ++i)
	{
		for (int j = i; j < r; ++j)
		{
			K(j,i) = pow(X(j) - X(i),2.0);
			K(i,j) = K(j,i);
		}
	}
	
}

double log_likelihood_micro2(const curve &Dat, const pq_point &theta, int k)
{

	int r = Dat.X.size();
	mat cov(r, r), SquareX(r, r);

	// caculate cov

	// vec colones = ones<vec>(r);
	// rowvec rowones = ones<rowvec>(r);
	mat EYE = eye<mat>(r, r);
	xixj(SquareX,Dat.X);
	// testfile<<SquareX;
	// square(Dat.X * rowones - colones * trans(Dat.X));

	cov = theta.v(k) * exp(-SquareX * theta.w(k) / 2.0) + EYE * theta.sigma2(k);
	// testfile << cov << endl;
	vec tmp(1);
	// testfile<<(cov * solve(cov, Dat.Y))<<endl;
	tmp = -trans(Dat.Y) * solve(cov, Dat.Y) / 2.0 - (double)r / 2.0 * log(2.0 * PI) - 1.0 / 2.0 * log(fabs(det(cov)));

	return tmp(0);
}
// more than one operator "/" matches these operands:
//  -- function template "arma::enable_if2<arma::is_arma_type<T1>::value,
//  const arma::eOp<T1, arma::eop_scalar_div_post>>::result arma::operator/(const T1 &X, T1::elem_type k)" 
// -- function template "arma::enable_if2<arma::is_arma_type<T1>::value, 
// const arma::eOp<T1, arma::eop_scalar_div_post>>::result operator/(const T1 &X, T1::elem_type k)" 
// -- operand types are: const arma::eOp<arma::eOp<arma::mat, arma::eop_neg>, arma::eop_scalar_times> / double

double log_likelihood2(const curve Data[], const pq_point &m)
{
	int t, k;
	int K = m.w.size();
	double norm;
	vec logP(K);
	norm = 0.0;
	pq_point theta;
	for (t = 0; t < Curve_num; t++)
	{
		logP = zeros<vec>(r);
		for (k = 0; k < K; k++)
		{
			// write_file(string s, ofstream &myfile)
			theta = m;
			// testfile << Data[t].Y << endl;
			logP(k) = log_likelihood_micro2(Data[t], m, k);
		}
		logP += log(m.pi);
		norm += logsumexp(logP);
		// testfile<<logP<<endl;
	}
	return norm;
}

double logsumexp(vec & logP){
	return log(sum(exp(logP - logP.max()))) + logP.max();
}


/************************************************************************/
/* Initializes the parameters of the model (and returns the scaling	*/
/* needed  to obtain a non -Inf log-likelihood, this should normaly be	*/
/* 1)									*/
/************************************************************************/

void draw_initial_model(const curve Data[], pq_point &m, double *logl)
{

	/* Start by drawing from prior					*/

	m.draw_from_prior(g);

	*logl = log_likelihood2(Data, m);
}

/************************************************************************/
/* Proposes a split move and returns log-likelihood			*/
/* put new classes in current and last position*/
/************************************************************************/
double prop_split(const curve Data[], pq_point &m, int k, int *k1, int *k2)
{
	double u[4];
	double logl;
	pq_point Elem_mix1(1), Elem_mix2(1);
	int K = m.w.size();

	if (K >= Curve_num)
	{

		exit(1);
	}
	if ((k < 0) || (k >= K))
	{
		exit(1);
	}

	u[0] = Beta(2, 2, g);
	u[1] = Beta(2, 2, g);
	u[2] = Beta(2, 2, g);
	u[3] = Beta(2, 2, g);

	/*_____________________________reallocate theta___________________*/
	double pi_k = m.pi(k);
	Elem_mix1.pi(0) = u[0] * m.pi(k);
	Elem_mix2.pi(0) = (1.0 - u[0]) * m.pi(k);

	Elem_mix1.sigma2(0) = u[1] * m.sigma2(k) * pi_k / Elem_mix1.pi(0);
	Elem_mix2.sigma2(0) = (1.0 - u[1]) * m.sigma2(k) * pi_k / Elem_mix2.pi(0);

	Elem_mix1.v(0) = u[2] * m.v(k) * pi_k / Elem_mix1.pi(0);
	Elem_mix2.v(0) = (1.0 - u[2]) * m.v(k) * pi_k / Elem_mix2.pi(0);

	Elem_mix1.w(0) = (1.0 - u[3]) / u[2] * m.w(k);
	Elem_mix2.w(0) = u[3] / (1.0 - u[2]) * m.w(k);

	/*____________________________insert w[k]________________________*/
	// m.deleteP_seq(k);
	double tmp = m.w(k);
	*k1 = m.locate_seq(Elem_mix1.w(0));
	m.insertPre_seq(*k1, Elem_mix1);
	*k2 = m.locate_seq(Elem_mix2.w(0));
	m.insertPre_seq(*k2, Elem_mix2);
	if ((*k1) < k)
		k++;
	if ((*k2) < k)
		k++;
	m.deleteP_seq(k);

	if ((*k1) > k)
		*k1 = (*k1) - 1;
	if ((*k2) > k)
		*k2 = (*k2) - 1;

	m.pi = normalise((m.pi), 1);

	logl = log_likelihood2(Data, m);
	return logl;
}

/************************************************************************/
/* Proposes a merge move and returns log-likelihood for new parameters  */
/* (puts the merged class in last position)				*/
/************************************************************************/
double prop_merge(const curve Data[], pq_point &m_old, pq_point &m_new, int *k, int k1, int k2)
{
	int i_old, i_new;
	double logl;
	pq_point Elem_mix(1);
	if (m_old.w.size() <= 1)
	{

		exit(1);
	}
	if ((k1 < 0) || (k2 < 0) || (k1 >= m_old.w.size()) || (k2 >= m_old.w.size()) || (k1 == k2))
	{

		exit(1);
	}

	/* Merge classes k1 and k2 in the last one of the new model		*/
	Elem_mix.pi(0) = m_old.pi(k1) + m_old.pi(k2);

	Elem_mix.sigma2(0) = (m_old.pi(k1) * m_old.sigma2(k1) +
						  m_old.pi(k2) * m_old.sigma2(k2)) /
						 Elem_mix.pi(0);

	Elem_mix.v(0) = (m_old.pi(k1) * m_old.v(k1) +
					 m_old.pi(k2) * m_old.v(k2)) /
					Elem_mix.pi(0);

	Elem_mix.w(0) = (m_old.pi(k1) * m_old.v(k1) * m_old.w(k1) +
					 m_old.pi(k2) * m_old.v(k2) * m_old.w(k2)) /
					Elem_mix.pi(0) / Elem_mix.v(0);

	/*_______________________merge_________________________*/
	m_new = m_old;

	*k = m_new.locate_seq(Elem_mix.w(0));
	m_new.insertPre_seq(*k, Elem_mix);

	if ((*k) < k1)
		k1++;
	if ((*k) < k2)
		k2++;
	m_new.deleteP_seq(k1);
	m_new.deleteP_seq(k2);

	if ((*k) > k1)
		*k = (*k) - 1;
	if ((*k) > k2)
		*k = (*k) - 1;

	m_new.pi = normalise(m_new.pi, 1);

	/* Compute new log-likelihood						*/
	logl = log_likelihood2(Data, m_new);
	return logl;
}

/************************************************************************/
/* Computes the part of the split ratio */
/************************************************************************/
double compute_log_split_ratio(pq_point &m, pq_point &m_split, int k, int k_split1, int k_split2)
{
	double lograt = 0.0;
	double u[4];
	int K = m.w.size();

	/*_________________________ Weights,pi__________________________*/
	u[0] = m_split.pi(k_split1) / m.pi(k);

	lograt += log((double)K) - log(6.0) - log(u[0] * (1 - u[0])) + log(m.pi(k));

	/*________________________ sigmav2s_____________________________*/
	u[1] = m_split.sigma2(k_split1) * m_split.pi(k_split1) / m.sigma2(k) / m.pi(k);
	boost::math::inverse_gamma_distribution<> IG_dist_s(10, 0.01);
	lograt += log(pdf(IG_dist_s, m_split.sigma2(k_split1))) + log(pdf(IG_dist_s, m_split.sigma2(k_split2))) - log(pdf(IG_dist_s, m.sigma2(k)));

	-log(6.0) - log(u[1] * (1 - u[1])) + log(m.sigma2(k)) + 2.0 * log(m.pi(k)) - log(m_split.pi(k_split1) * m_split.pi(k_split2));

	/*________________________v0s________________________________*/
	u[2] = m_split.v(k_split1) * m_split.pi(k_split1) / m.v(k) / m.pi(k);
	boost::math::inverse_gamma_distribution<> IG_dist_v(0.5, 0.5);
	lograt += log(pdf(IG_dist_v, m_split.v(k_split1))) + log(pdf(IG_dist_v, m_split.v(k_split2))) - log(pdf(IG_dist_v, m.v(k)));

	-log(6.0) - log(u[2] * (1 - u[2])) + log(m.v(k)) + 2.0 * log(m.pi(k)) - log(m_split.pi(k_split1) * m_split.pi(k_split2));

	/*________________________w________________________________*/
	u[3] = m_split.w(k_split2) / m.w(k) * (1 - u[2]);
	boost::math::inverse_gamma_distribution<> IG_dist_w(0.5, 0.5);
	lograt += log(pdf(IG_dist_w, m_split.w(k_split1))) + log(pdf(IG_dist_w, m_split.w(k_split2))) - log(pdf(IG_dist_w, m.w(k)));

	-log(6.0) - log(u[3] * (1 - u[3])) + log(m.w(k)) - log(u[2] * (1 - u[2]));

	return lograt;
}

double calc_secondary_moment(const curve &Dat, const pq_point &m, const pq_point &m_new, int k, int k1, int k2)
{
	double secondary_moment = 0.0;
	// evaluate the covariance for the mth batch/curve
	// there are totally Nm points on this batch/curve
	int r = Nm;
	mat cov(r, r), cov1(r, r), cov2(r, r), DiffCov(r, r), SquareX(r, r);

	// vec colones = ones<vec>(r);
	// rowvec rowones = ones<rowvec>(r);
	mat EYE = eye<mat>(r, r);

	xixj(SquareX,Dat.X);
	//  square(Dat.X * rowones - colones * trans(Dat.X));

	cov = m.v(k) * exp(-SquareX * m.w(k) / 2.0);
	cov = cov + EYE * m.sigma2(k);
	cov1 = m_new.v(k1) * exp(-SquareX * m_new.w(k1) / 2.0) + EYE * m_new.sigma2(k1);
	cov2 = m_new.v(k2) * exp(-SquareX * m_new.w(k2) / 2.0) + EYE * m_new.sigma2(k2);

	DiffCov = m.pi(k) * cov - m_new.pi(k1) * cov1 - m_new.pi(k2) * cov2;

	return secondary_moment = norm(DiffCov, "fro");
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
