
#include <string.h>
#include "mix_lib.h"
#include <iostream>
#include <fstream>
#include "armadillo"
#include <ctime>

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
double w_eps, v0_eps, sigmav2_eps;


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

								w_eps = 5e-6;
								v0_eps = 0.5;
								sigmav2_eps = 0.005;
								printf("Data has %d curves. Running %d x %d iterations of the sampler...\n", Curve_num, NOut, SubSamp);
}

int main(int argc, char **argv)
{
								curve Data[MaxM];
								pq_point theta(3);
								// double ran;
								double logl;
								read_parameters(argc, argv, Data);
								// cout<<Data[0].X;
								mat K(Nm,Nm);
								time_t t_start, t_end;
								// // ofstream TimeFP("time.txt");
								double DiffTime;
								t_start = time(NULL);
								for(int i=0; i<10000; i++)
																xixj(K,Data[0].X);
								t_end = time(NULL);
								DiffTime = difftime(t_end, t_start);
								cout<<K;
								cout<<DiffTime<<endl;

}
