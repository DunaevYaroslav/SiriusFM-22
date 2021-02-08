#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <thread>
namespace SiriusFM{
    enum class CcyE
    {
        UNDEFINED = -1,
        USD = 0,
        EUR = 1,
        GBP = 2,
        CHF = 3,
        RUB = 4,
        N = 5
    };
   inline char const*CcyE2Str(CcyE a_ccy)
    {
        switch(a_ccy)
            {
                case CcyE::USD:return "USD";
                case CcyE::RUB:return "RUB";
                case CcyE::GBP:return "GBP";
                case CcyE::CHF:return "CHF";
                case CcyE::EUR:return "EUR";
                default: throw std::invalid_argument("Bad arg");
            }
    }
    inline CcyE Str2CcyE(char const*a_str)
    {
        if(strcmp(a_str, "USD")==0){
            return CcyE::USD;
        }
            if(strcmp(a_str, "RUB")==0){
                 return CcyE::RUB;
            }
            if(strcmp(a_str, "EUR")==0){
                 return CcyE::EUR;
        }
                if(strcmp(a_str, "GBP")==0){
            return CcyE::GBP;
        }
                if(strcmp(a_str, "CHF")==0){
            return CcyE::CHF;
        }
        else
        {
            throw std::invalid_argument("Bad arg");
        }
        
    }
    enum class IRmodeE
    {
        Const = 0,
        FwdCurve = 1,
        Stoch = 2
    };
    template<IRmodeE IRM>
    class IRProvider;
}