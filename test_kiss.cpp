
#include <string.h>
// #include "nuts.h"
#include <iostream>
#include <fstream>
#include "armadillo"
#include "mix_lib.h"
// #define ULONG_NORM 18446744073709551616.0
/*Definitions                */
#define NArgs 17
#define StrLen 40
#define ArgExt ".arg"
#define ResExt ".res"
#define StatExt ".sts"

//
//using namespace std;
//using namespace arma;
/*Global variables                            */
/*Seeds for the kiss generator (see alea)     */

/*Files                                      */
Generator g;
char DataFile[StrLen], StatFile[StrLen];
ofstream parameterFP, zFP, StatFP, testfile;
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
int Curve_num;
int Nm;
using namespace std;
using namespace arma;
struct STATS
{
	string split_or_merge, acc_or_rej, simu_acc_or_simu_rej, delete_empty_component;
	double sm_prob, simu_prob;
	void initialize()
	{
		split_or_merge = "NULL";
		acc_or_rej = "NULL";
		simu_acc_or_simu_rej = "NULL";
		delete_empty_component = "NULL";
		sm_prob = 0.0;
		simu_prob = 0.0;
	}
	void print()
	{
		cout << split_or_merge
			 << " " << acc_or_rej << " " << sm_prob << " " << simu_acc_or_simu_rej << " " << simu_prob << " " << delete_empty_component << endl;
	}
	void write_file(ofstream &myfile)
	{
		myfile << split_or_merge
			   << " " << acc_or_rej << " " << sm_prob << " " << simu_acc_or_simu_rej << " " << simu_prob << " " << delete_empty_component << endl;
	}
};

void read_parameters(int argc, char **argv, curve Data[])
{
        testfile.open("generator.res");
        char argfile[StrLen];

        strcpy(StatFile, argv[1]);
        strcpy(DataFile, argv[1]);
        strcpy(argfile, argv[1]);

        strcat(argfile, ArgExt);
        ifstream argfp(argfile);
        argfp >> g.i >> g.j >> g.k >> NOut >> SubSamp >> Kmax >> PFixed >> PBirth >> PDeath >> PSplit;
        argfp.close();
        /* Compute frequently used quantities					*/
        NIt = NOut * SubSamp;
        PMerge = 1.0 - (PFixed + PBirth + PDeath + PSplit);
        PFixed_or_BD = PFixed + PBirth + PDeath;

        /* Read data								*/
        strcat(DataFile, ".dat");
        mat AAA;
        AAA.load(DataFile);
        Curve_num = AAA.n_rows / 2;
        Nm = AAA.n_cols;
        for (int m = 0; m < Curve_num; m++)
        {
                Data[m].X = AAA.row(2 * m).t();
                Data[m].Y = AAA.row(2 * m + 1).t();
        }

        // w_eps = 5e-6;
        // v0_eps = 0.5;
        // sigmav2_eps = 0.005;
        printf("Data has %d curves. Running %d x %d iterations of the sampler...\n", Curve_num, NOut, SubSamp);
}

/************************************************************************/
/* Write iterations on disk						*/
/************************************************************************/
void write_data(pq_point &theta, double logl, STATS &stats, vec &z, int in)
{

        
                /* Initial value, we need to open the files				*/
                /* Open statfile							*/
                strcat(StatFile, StatExt);
                StatFP.open(StatFile);

                parameterFP.open("parameter.res");
                zFP.open("z.res");
        

        
                StatFP << in << " " << theta.w.size() << " " << logl << " ";
                stats.write_file(StatFP);
              
                theta.write_file(parameterFP);
                zFP << z.t();
        
                StatFP.close();
                parameterFP.close();
                zFP.close();
        
}


int main(int argc, char **argv)
{
	curve Data[MaxM];
	pq_point theta(3);
	double ran;
	double logl;
	int in;
        vec logP(2);

        logP<<1<<1<<endr;
        double c=logsumexp(logP);
        cout<<c<<endl;
	vec z;
	STATS stats;
	z << 0 << 0 << 0 << 1 << 1 << 1 << 2 << 2 << 2 << endr;
	stats.initialize();
	read_parameters(argc, argv, Data);
	draw_initial_model(Data, theta, &logl);
	log_likelihood2(Data, theta);

        theta.print("theta:");
        write_data(theta, logl, stats, z, 1);

       
}
