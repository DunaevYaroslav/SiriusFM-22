#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <thread>
namespace SiriusFM{
    class DiffusionOU
    {
        private:
            double const const_sigma;
            double const teta;
            double const kappa;
        public:
            double mu(double St, double t) const{
                return kappa*(teta-St);
            }
            double sigma(double St, double t) const{
                if(const_sigma<=0){return 0;}
                if(kappa<0){return 0;}
                if(teta<=0){return 0;}
                return const_sigma*sqrt(St);
            }
        DiffusionOU(double const_s, double const_teta, double const_kappa)
        :const_sigma(const_s),
         kappa(const_kappa),
         teta(const_teta){
        if(const_sigma<=0 || kappa < 0 || teta <= 0){throw std::invalid_argument("Bad arguments");}
         }
    };
}