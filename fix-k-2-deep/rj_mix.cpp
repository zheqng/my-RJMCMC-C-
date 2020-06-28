
#include <string.h>
#include "nuts.h"
#include <iostream>
#include <fstream>
// #include "armadillo"
#include <ctime>



/*Definitions                */
#define NArgs 17
#define StrLen 40
#define ArgExt  ".arg"
#define ResExt  ".res"
#define StatExt  ".sts"
#define MAXBUFSIZE  ((int) 1e6)


//
//using namespace std;
//using namespace arma;
/*Global variables                            */
/*Seeds for the kiss generator (see alea)     */

/*Files                                      */
Generator g;
char DataFile[StrLen], StatFile[StrLen];
ofstream parameterFP,zFP, StatFP,testfile;
/*Prior hyperparameters                    */
int Kmax;
/*Sampler settings                        */
int NOut, SubSamp, NIt;

/*Fixed k move                             */
double PFixed;
/*Birth and death                          */
double PBirth, PDeath, PFixed_or_BD;
/*Split and merge                          */
double PSplit, PMerge;
int Curve_num, THREAD_NUM;
int Nm;
double w_eps,v0_eps,sigmav2_eps;



/* Function prototypes							*/
struct STATS {
        string split_or_merge,acc_or_rej,simu_acc_or_simu_rej,delete_empty_component;
        double sm_prob,simu_prob;
        void initialize(){
                split_or_merge = "NULL";
                acc_or_rej = "NULL";
                simu_acc_or_simu_rej = "NULL";
                delete_empty_component = "NULL";
                sm_prob = 0.0;
                simu_prob = 0.0;
        }
        void print(int flag){
                if(flag==1) //print split or merge information
                        cout <<  split_or_merge
                             <<" "<<acc_or_rej<<endl;
                // "; simulated annealing "<<
                //    simu_acc_or_simu_rej<<endl;
                else         //print delete empty component information
                //<<" simulated annealing "<<
                //  simu_acc_or_simu_rej<<"; "
                {
                        if(delete_empty_component=="delete")
                                cout<<delete_empty_component<<" empty components"<<endl;
                }


        }
        void write_file(ofstream & myfile){
                myfile<<  split_or_merge
                      <<" "<<acc_or_rej<<" "<<sm_prob<<" "<<
                        simu_acc_or_simu_rej<<" "<<simu_prob<<" "<<delete_empty_component<<endl;
        }
};
void read_parameters(int argc, char** argv, curve Data[] );
MatrixXd readMatrix(const char *filename);
void write_data(pq_point & theta, double logl,STATS & stats,VectorXi & z,int in);
void RJMH_birth_or_death(curve Data[], pq_point & m,double* logl,VectorXi & z,STATS & stats);
void RJMH_split_or_merge(curve Data[], pq_point & m, double* logl,STATS & stats);
double Simulated_Annealing( double *logl, double *logl_old,pq_point & theta, pq_point & theta_old,int in,STATS & stats);
void  sample_nuts_cpp(const curve Data[], pq_point & current_theta,
                      VectorXi & z);
void write_theta(pq_point & theta, double logl,VectorXi & z,int in);

int main(int argc, char** argv) {
        curve Data[MaxM];
        pq_point theta(2);
        double ran;
        double logl;
        int in;
        VectorXi z;
        STATS stats;

        stats.initialize();

        read_parameters(argc, argv, Data);

        draw_initial_model(Data, theta, &logl);
        theta.print("theta:");
        Gibbs_Sampling_z(Data, theta,z);
        write_data(theta, logl, stats,z,0);



        /* Main loop								*/
        pq_point theta_old;
        double logl_old,ratio= 1.0;
        for (in=1; in<=NIt; in++) {
                cout<<in<<"th iter"<<endl;
                /* Fixed k move							*/
                /*_____________update pi________________________________*/
                // if( kiss(g)<0.0) {
                //         RJMH_split_or_merge(Data, theta, &logl, stats);
                //         // // if(stats.simu_acc_or_simu_rej=="accept")
                //         // {  stats.print(1);
                //         //
                //         //    cout<<theta.v<<endl;}
                // }
                // else{
                /*************sample pi************************/
                cout<<"gibbs sample z:"<<endl;
                Gibbs_Sampling_z(Data,theta,z);
                cout<<"gibbs sample pi"<<endl;
                Gibbs_Sampling_pi(theta,z);
                cout<<"nuts"<<endl;
                /*___________________sample theta__________________*/
                sample_nuts_cpp(Data,theta,z);
                theta.print("theta:");

                // cout<<"after nuts, v:"<<endl;
                // cout<<theta.v<<endl;
                // cout<<"gibbs sample pi:"<<endl;
                Gibbs_Sampling_pi(theta,z);
                // cout<<"gibbs sample z:"<<endl;
                Gibbs_Sampling_z(Data,theta,z);
                /*_________________change pi and K____________________________*/
                // cout<<"delete empty component:"<<endl;
                // RJMH_birth_or_death(Data, theta, &logl,z,stats);
                // theta.print("theta:");

                //split_or_merge,acc_or_rej,simu_acc_or_simu_rej,delete_empty_component;
                // if(stats.simu_acc_or_simu_rej=="accept")
                // {stats.print(0);
                //  // cout<<"after "<<stats.split_or_merge<<stats.acc_or_rej<<", v:"<<endl;
                //  cout<<theta.v<<endl;}
                // }

                // if(in>1) {
                //         // cout<<2<<endl;
                //         ratio = Simulated_Annealing(&logl,&logl_old,theta,theta_old,in,stats);
                //         theta.print("theta:");
                //
                // }

                theta_old = theta;
                logl_old = logl;
                /* Note: The test below is true for in = 0 and in = NIt		*/
                //if ((div(in, SubSamp)).rem == 0){
                // Gibbs_Sampling_z(Data,theta,z);
                // write_data(theta, logl, stats,z, in);


        }

        testfile.close();
        return(0);
}


MatrixXd readMatrix(const char *filename)
{
        int cols = 0, rows = 0;
        double buff[MAXBUFSIZE];

        // Read numbers from file into buffer.
        ifstream infile;
        infile.open(filename);
        while (!infile.eof())
        {
                string line;
                getline(infile, line);

                int temp_cols = 0;
                stringstream stream(line);
                while(!stream.eof())
                        stream >> buff[cols*rows+temp_cols++];

                if (temp_cols == 0)
                        continue;

                if (cols == 0)
                        cols = temp_cols;

                rows++;
        }

        infile.close();
        // cout<<c/ols<<endl;
        rows--;

        // Populate matrix with numbers.
        MatrixXd result(rows,cols);
        for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                        result(i,j) = buff[ cols*i+j ];

        return result;
};
/************************************************************************/
/* Read parameters (including seeds for the random generator) and data  */
/************************************************************************/
void read_parameters(int argc, char** argv, curve Data[] ) {
        testfile.open("generator.res");
        char argfile[StrLen];

        strcpy(StatFile, argv[1]);
        strcpy(DataFile, argv[1]);
        strcpy(argfile, argv[1]);

        strcat(argfile, ArgExt);
        ifstream argfp(argfile);
        argfp>>g.i>>g.j>>g.k>>NOut>>THREAD_NUM>>SubSamp>>Kmax>>PFixed>>PBirth>>PDeath>>PSplit;
        argfp.close();
        /* Compute frequently used quantities					*/
        NIt = NOut*SubSamp;
        PMerge = 1.0 - (PFixed + PBirth + PDeath + PSplit);
        PFixed_or_BD = PFixed + PBirth + PDeath;

        /* Read data								*/
        strcat(DataFile,".dat");
        MatrixXd AAA =  readMatrix(DataFile);


        Curve_num = AAA.rows()/2;
        Nm = AAA.cols();

        for(int m=0; m<Curve_num; m++) {
                Data[m].X = AAA.row(2 * m).segment(0, Nm - 1).transpose();
                Data[m].Y = AAA.row(2 * m + 1).segment(0, Nm - 1).transpose();
        }
        // cout<<Data[0].X<<endl;
        // cout<<Data[0].Y<<endl;
        w_eps = 5e-6;
        v0_eps = 0.5;
        sigmav2_eps = 0.005;
        printf( "Data has %d curves. Running %d x %d iterations of the sampler...\n", Curve_num, NOut, SubSamp);
        cout<<THREAD_NUM<<endl;
}



/************************************************************************/
/* Write iterations on disk						*/
/************************************************************************/
void write_data(pq_point & theta, double logl, STATS & stats,VectorXi & z,int in) {

        if (in == 0) {
                /* Initial value, we need to open the files				*/
                /* Open statfile							*/
                strcat(StatFile, StatExt);
                StatFP.open(StatFile);

                parameterFP.open("parameter.res");
                zFP.open("z.res");
        }

        else{
                StatFP<<in<<" "<<theta.w.size()<<" "<<logl<<" ";
                stats.write_file(StatFP);
                // <<stats.split_or_merge
                // <<" "<<stats.acc_or_rej<<" "<<stats.sm_prob<<" "<<
                // stats.simu_acc_or_simu_rej<<" "<<stats.simu_prob<<endl;
                // ResFP<<in<<endl;
                theta.write_file(parameterFP);
                zFP<<z.transpose()<<endl;
        }
        if(in == NOut) {
                StatFP.close();
                parameterFP.close();
                zFP.close();
        }
}


void RJMH_birth_or_death(curve Data[], pq_point & m,double* logl,VectorXi & z,STATS & stats) {
        int K=m.w.size();
        VectorXi label_num=VectorXi::Zero(K);
        int k;

        for(int i = 0; i< Curve_num; i++) {k = z(i); label_num(k)+=1;}
        // cout<<"number of class:"<<label_num<<endl;
        VectorXi empty_component(K);
        int j=0, empty_component_num=0;
        for(int k=0; k<K; k++) {
                if(label_num(k)==0)
                {empty_component(j)=k; j++;}
        }
        // cout<<"empty component is:"<<empty_component<<endl;
        empty_component_num = j;
        // cout<<"the number of empty component is:"<<empty_component_num<<endl;
        // cout<<"before delete:"<<m.v<<endl;
        for(int i=0; i<empty_component_num; i++) {
                m.deleteP_seq(empty_component(i)-i);
        }
        // cout<<"after delete:"<<m.v<<endl;
        m.pi = m.pi/m.pi.sum();

        if(empty_component_num>0) {
                stats.delete_empty_component = "delete";
        }
        else stats.delete_empty_component = "reserve";
        *logl = log_likelihood2(Data,m);
}

/************************************************************************/
/* Implements the RJ MH move based on split/merge proposal		*/
/************************************************************************/
void RJMH_split_or_merge(curve Data[], pq_point & m, double* logl, STATS & stats) {
        double prop_ratio, add_logratio,ratio;
        double logl_new;
        pq_point m_new;
        int k, k1, k2,K=m.w.size();
        double split, accept;
        VectorXi z_new;
        double secondary_moment;
        /* Take care of case 1, 2, M-1 and M components (this could be done	*/
        /* more simply)							*/
        if (K == 1) {
                split = 1;

                prop_ratio = PMerge/(1.0-PFixed_or_BD);
        }
        else if (K == Kmax) {
                split = 0;
                prop_ratio = PSplit/(1.0-PFixed_or_BD);
        }
        else {
                if (kiss(g) < PSplit/(1.0-PFixed_or_BD)) {
                        /* Split								*/
                        split = 1;
                        if (K == (Kmax-1))
                                prop_ratio = (1.0-PFixed_or_BD)/PSplit;
                        else
                                prop_ratio = PMerge/PSplit;
                }
                else {
                        /* Merge								*/
                        split = 0;
                        if (K == 2)
                                prop_ratio = (1.0-PFixed_or_BD)/PMerge;
                        else
                                prop_ratio = PSplit/PMerge;
                }
        }


        if(split) {
                cout<<"split:"<<endl;
                stats.split_or_merge = "split";
                m_new = m;
                k = (int)floor((double)K *kiss(g));
                if(k==K) k--;
                /* Proposes a split move and returns log-likelihood			*/

                logl_new = prop_split(Data, m_new, k,&k1,&k2);
                secondary_moment = calc_secondary_moment((Data[0]),m,m_new,k,k1,k2);
                add_logratio = compute_log_split_ratio(m, m_new, k, k1,k2);
        }
        else {
                cout<<"merge"<<endl;
                /* Draw two distinct indices using modulo m->k addition		*/
                stats.split_or_merge = "merge";
                k1 = (int)floor((double)(K-1) * kiss(g));
                if(k1==K-1) k1--;
                k2 = 1+k1;
                // cout<<"k1:"<<k1<<" k2:"<<k2<<endl;
                logl_new = prop_merge(Data, m, m_new, &k,k1, k2);
                // cout<<"k:"<<k<<"k1:"<<k1<<" k2:"<<k2<<endl;
                // cout<<3<<endl;
                add_logratio =-compute_log_split_ratio(m_new, m, k, k1, k2);
                // cout<<4<<endl;
                // secondary_moment = calc_secondary_moment((Data[0]),m_new,m,k,k1,k2);
        }

        double x;
        x = kiss(g);
        ratio = exp((logl_new-(*logl))+add_logratio);
        ratio *= prop_ratio;
        if(ratio>1) ratio =1.0;
        if(ratio>0) {}
        else
                ratio = 0.0;
        // if (secondary_moment>0.01) ratio= -1.0;
        /* Accept/reject							*/
        if (x<ratio) {
                accept = 1;
                /* Modify the parameters and log likelihood				*/
                m=m_new;
                *logl = logl_new;
                stats.acc_or_rej = "accept";
                // cout<<"secondary_moment "<<secondary_moment<<endl;
        }
        else{
                accept = 0;
                stats.acc_or_rej= "reject";
        }

        stats.sm_prob = ratio;

}


double Simulated_Annealing( double *logl, double *logl_old,pq_point & theta, pq_point & theta_old,int in,STATS & stats){
        int K = theta.w.size();
        int old_K = theta_old.w.size();

        double add_logratio = compute_log_prior_ratio(theta,theta_old);

        double criation = ( 4.0 *K)*log((double)(Curve_num*Nm));
        double criation_old = ( 4.0 *old_K)*log((double)(Curve_num*Nm));

        double Ta = 1.0, Tf = 1e-5;
        double Tb = (Ta - Tf)/(double)NIt;
        double x = kiss(g);
        double ratio = (((*logl)-(*logl_old))+add_logratio -criation/2.0+criation_old/2.0);

        // testfile<<in<<"th iter: Curve num="<<Curve_num<<" Nm="<<Nm<<" criation="<<criation
        // <<" criation_old="<<criation_old<<endl;
        // testfile<<"add_logratio="<<add_logratio<<" logl="<<*logl<<" logl_old="<<*logl_old<<endl;
        double temp = Ta-Tb*in;
        //tempureture from 1 to 1e-5
        ratio = ratio;
        ratio = exp(ratio);
        if(ratio>1.0) ratio =1.0;
        if(ratio>0.0) {}
        else
                ratio = 0.0;
        stats.simu_prob = ratio;
        if (x<ratio) {
                stats.simu_acc_or_simu_rej = "accept";
        }
        else{
                theta = theta_old;
                *logl = *logl_old;
                stats.simu_acc_or_simu_rej = "reject";
        }
        return ratio;
}



void sample_nuts_cpp(const curve Data[], pq_point &current_theta,
                     VectorXi &z)
{


        int iter = 1000;
        int MAXDEPTH = 8;
        double lambda = 1.0;
        int KK = current_theta.w.size();

        int I_adapt = iter/2.0, t0 = iter;
        double delta = 0.65, kappa = 0.75, gamma = 0.05;
        double H_bar[100], log_eps_bar[100], log_eps[100];
        double epsilon = 1e-6;
        log_eps[0] = log(epsilon);
        H_bar[0] = 0;
        double mu = log(10.0 * epsilon);
        double logl;
        time_t t_start,t_end;
        ofstream TimeFP("time.txt");
        double DiffTime;
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
        t_start = time(NULL);
        t_end = time(NULL);
        DiffTime = difftime(t_end,t_start);
        TimeFP<<"0 "<<DiffTime<<endl;
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
                        cout<<"build a new tree, depth="<<depth_<<endl;
                        // testfile << "p_sharp_plus:" << endl;
                        // testfile<< p_sharp_plus.w.transpose() << endl;
                        // testfile << p_sharp_plus.v.transpose() <<endl;
                        // testfile << p_sharp_plus.sigma2.transpose()<<endl;
                        //
                        // testfile << "p_sharp_minus:" << endl;
                        // testfile<< p_sharp_minus.w.transpose() << endl;
                        // testfile << p_sharp_minus.v.transpose() <<endl;
                        // testfile << p_sharp_minus.sigma2.transpose()<<endl;



                        // Build a new subtree in a random direction
                        util.sign = 2 * (kiss(g) < 0.5) - 1;
                        int n_valid_subtree = 0;
                        if (util.sign == 1)
                        {
                                phi=phi_plus;
                                theta=theta_plus;
                                n_valid_subtree = BuildTree(Data, phi, theta, phi_propose, theta_propose,
                                                            p_sharp_dummy, p_sharp_plus, rho_subtree, util, depth_, epsilon, lambda, z);
                                phi_plus=phi_propose;
                                theta_plus=theta_propose;
                        }
                        else
                        {
                                phi = phi_minus;
                                theta = theta_minus;
                                n_valid_subtree = BuildTree(Data, phi, theta, phi_propose, theta_propose,
                                                            p_sharp_minus, p_sharp_dummy, rho_subtree, util, depth_, epsilon, lambda, z);
                                phi_minus=phi_propose;
                                theta_minus=theta_propose;
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
                Gibbs_Sampling_z(Data,current_theta,z);
                logl =  log_likelihood2( Data, current_theta);
                write_theta(current_theta, logl, z, i);
                t_end = time(NULL);
                DiffTime = difftime(t_end,t_start);
                TimeFP<<i<<" "<<DiffTime<<" ";
                cout<<i<<" "<<DiffTime<<" ";
                // cout<<"Difftime:"<<DiffTime<<"log likelihood:"<<logl<<endl;
                if (i%10 == 0) TimeFP<<endl;
        } // end for
        TimeFP.close();
        // return samples;
}



void write_theta(pq_point & theta, double logl,VectorXi & z,int in) {

        if (in == 0) {
                /* Initial value, we need to open the files				*/
                /* Open statfile							*/
                // strcat(StatFile, StatExt);
                // StatFP.open(StatFile);

                parameterFP.open("parameter.res");
                zFP.open("z.res");
        }

        else{
                // <<stats.split_or_merge
                // <<" "<<stats.acc_or_rej<<" "<<stats.sm_prob<<" "<<
                // stats.simu_acc_or_simu_rej<<" "<<stats.simu_prob<<endl;
                // ResFP<<in<<endl;
                theta.write_file(parameterFP);
                zFP<<z.transpose()<<endl;
        }
        if(in == 1000) {

                parameterFP.close();
                zFP.close();
        }
}
