#include <iostream>
#include <string>
#include <omp.h>
#include "mix_lib.h"

using namespace std;
void test( ) {
        int a=0;
//      cout<<"omp_get_max_threads()"<<omp_get_max_threads()<<" omp_get_num_threads() "<<omp_get_num_threads()<<endl;
        omp_set_nested(1);
  #pragma omp parallel for
        for(int i=0; i<55; i++)
        {
                // cout<<"The outer loop:thread_num"<<omp_get_thread_num()<<endl;
  #pragma omp critical
                {
          #pragma omp parallel for
                        for(int j=0; j<55; j++)
                        {
                                a +=i+j;
                                cout<<"the inner loop: thread_num"<<omp_get_thread_num()<<endl;
                        }
                }
        }

}
int main(){

        // test();


        return 0;
}
