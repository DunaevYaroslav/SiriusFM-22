#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
namespace SiriusFM{
    class DiffusionGBM
    {
        private:
            double const const_mu;
            double const const_sigma;
        public:
            double mu(double St, double t) const{
                return const_mu*St;
            }
            double sigma(double St, double t) const{
                if(const_sigma<=0){return 0;}
                return const_sigma*St;
            }
        DiffusionGBM(double const_m, double const_s)
        :const_mu(const_m),
         const_sigma(const_s)
        {
            if(const_sigma<=0){throw std::invalid_argument("Bad arguments");}
        }
    };
}