#include <array>
#include "utils.h"

namespace particleKernels
{

template<int dims>
struct twoBodyPotential
{
    using Real=double;
    static const int dimensions = dims;

    twoBodyPotential(std::array<int,2> rangeA,std::array<int,2> rangeB,std::array<Real,dims> lBox_ );
    
    void setPeriodic( std::array<Real,dims> lBox);
    
    template<class V_t>
    Real operator()(const V_t & V,const Real * positions, int i1 , int i2, int t0 , int t1,int N, int T) const;

    template<class V_t>
    void addForce(const V_t & V,const Real * positions,Real * forces, int i1 , int i2, int t0 , int t1,int N, int T) const;




    bool containedInSetA(int i1, int i2) const;
    bool containedInSetB(int i1, int i2) const;


    private:

    bool periodic=false;
    int iStartA, iStartB, iEndA,iEndB;
    std::array<Real,dimensions> _lBox;
    bool isTriangular=false;
};


#include "twoBodyPotential.hpp"


}


