#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <thread>
namespace SiriusFM{
    class DiffusionLipton
    {
        private:
            double const const_sigma0;
            double const const_sigma1;
            double const const_sigma2;
            double const const_mu;
        public:
            double mu(double St, double t) const{
                return const_mu*St;
            }
            double sigma(double St, double t) const{
                if(const_sigma0<=0 || const_sigma2<=0 || (pow(const_sigma1, 2.0) - 4*const_sigma0*const_sigma2>=0)){return 0;}
                return const_sigma0+const_sigma1*St+const_sigma2*pow(St, 2.0);
            }
        DiffusionLipton(double const_s0, double const_s1, double const_s2, double const_m)
        :const_sigma0(const_s0),
         const_sigma1(const_s1),
         const_sigma2(const_s2),
         const_mu(const_m){
        if(const_sigma0<=0 || const_sigma2<=0 || (pow(const_sigma1, 2.0) - 4*const_sigma0*const_sigma2>=0)){throw std::invalid_argument("Bad arguments");}
         }
    };
}
