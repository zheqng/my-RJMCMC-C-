#include <iostream>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>




#if __WORDSIZE == 64
#  define ULONG_NORM  18446744073709551616.0
#else
#  define ULONG_NORM  4294967296.0
#endif

using namespace boost::random;

using namespace std;

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
        // void write_file(string s,ofstream & myfile){
        //         myfile<<s<<i<<" "<<j<<" "<<k<<endl;
        // }
};

double kiss(Generator & g) {

        double x;

        (g.j) = (g.j) ^ ((g.j)<<17);
        (g.k) = ( (g.k) ^ ( (g.k) <<18)) & 0x7FFFFFFF;
        (g.i) = 69069UL*((g.i))+23606797UL;
        (g.j) ^= ((g.j)>>15);
        (g.k) ^= ((g.k)>>13);
        x = (double)((g.i)+(g.j)+(g.k)) / ULONG_NORM;
        // cout<<"kiss(g):"<<x<<" ";
        return(x);
}

int main(){
        Generator g0;
        g0.i = 113584957;
        g0.j = 487306701;
        g0.k = 812614905;



        // boost::random::mt19937 gen(time(0));
        // boost::uniform_01<boost::mt19937&> u01(gen);

        for(int i=0; i<1000; i++) {
                kiss(g0);
        }
        g0.print("g0:");
        cout<<kiss(g0)<<endl;
        cout<<ULONG_NORM<<endl;


}
