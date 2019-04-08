
#include <string.h>
// #include "nuts.h"
#include <iostream>
#include <fstream>
#include "armadillo"
#  define ULONG_NORM  18446744073709551616.0

using namespace std;

// #include <ctime>
struct Generator {
								unsigned long int i,j,k;
public:
								Generator& operator=(const Generator & g0){
																if(this == &g0)
																								return *this;

																i = g0.i;
																j = g0.j;
																k = g0.k;

																return *this;
								}
								void print(string s){
																cout<<s<<i<<" "<<j<<" "<<k<<endl;
								}
								void write_file(string s,ofstream & myfile){
																myfile<<s<<i<<" "<<j<<" "<<k<<endl;
								}
};


double kiss(Generator & g) {
								double x;

								(g.j) = (g.j) ^ ((g.j)<<17);
								(g.k) = ( (g.k) ^ ( (g.k) <<18)) & 0x7FFFFFFF;
								(g.i) = 69069UL*((g.i))+23606797UL;
								(g.j) ^= ((g.j)>>15);
								(g.k) ^= ((g.k)>>13);
								x = (double)((g.i)+(g.j)+(g.k)) / ULONG_NORM;
								return(x);
}

int main(){
								Generator g;
								g.i = 113584957;
								g.j = 487306701;
								g.k = 812614905;
								for(int i=0; i<100; i++) {
																cout<<kiss(g)<<" ";
																if((i%10)==0) cout<<endl;
								}
}
