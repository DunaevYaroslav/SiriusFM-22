#include <cmath>
#include <ctime>
#include<stdio.h>
#include"IRProvider.h"
#include"IRProviderConst.h"
#include"DiffusionOU.h"
#include"DiffusionLipton.h"
#include"DiffusionCEV.h"
#include"DiffusionCIR.h"
#include"DiffusionGBM.h"
namespace SiriusFM{
template <typename Diffusion1D, typename AProvider, typename BProvider, typename AssetClassA, typename AssetClassB>
class MCEngine1D
{
    private:
        long const m_MaxL;
        long const m_MaxP;
        double * const m_paths;
        long m_L;
        long m_P;
        double m_T;
        double m_t0;
        Diffusion1D const * m_diff;
        AProvider const * m_rateA;
        AProvider const * m_rateB;
        //AE m_ae;
        //BE m_be;
        bool m_isRN;
    public:
        MCEngine1D(long a_MaxL, long a_MaxP)
        :m_MaxL(a_MaxL),
        m_MaxP(a_MaxP),
        m_paths(new double[m_MaxL*m_MaxP]),
        m_L(0),
        m_P(0),
        m_T(nan),
        m_t0(nan),
        m_diff(nullptr),
        m_rateA(nullptr),
        m_rateB(nullptr),
        m_a(AssetClassA::UNDEFINED),
        m_b(AssetClassB::UNDEFINED),
        m_isRN(false){
            if(m_MaxL<=0 || m_MaxP<=0){
                throw std::invalid_argument("Bad arg");
            }
        }
        void Simulate(time_t a_t0, time_t a_T, int a_Tas_min, Diffusion1D const * a_diff, 
        AProvider const * a_rateA, BProvider const * a_rateB, AssetClassA a_A, AssetClassB a_B, bool a_isRN);
    };
}