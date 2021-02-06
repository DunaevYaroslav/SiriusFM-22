#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <thread>
#include"IRProvider.h"
namespace SiriusFM{
    template<>
    class IRProvider<IRmodeE::Const>
    {
        private:
        double m_IRS[int(CcyE::N)];
        public:
        IRProvider(char const * a_file);
        double r(CcyE a_ccy, double a_t) const
        {return m_IRS[int(a_ccy)];}
    };
}
