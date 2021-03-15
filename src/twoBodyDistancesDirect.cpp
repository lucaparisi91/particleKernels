#include "twoBodyDistancesDirect.h"
#include "utils.h"
#include <cassert>


using Real = double;
using namespace particleKernels;


template<int dims>
twoBodyPotential<dims>::twoBodyPotential(std::array<int,2> rangeA,std::array<int,2> rangeB, std::array<Real,dims> lBox_  ) 
{
    iStartA = rangeA[0];
    iStartB = rangeB[0];
    
    iEndA = rangeA[1];
    iEndB = rangeB[1];

    setPeriodic(lBox_);


    // check that the two ranges do not intersect
    assert(iStartA<=iEndA);
    assert(iStartB<=iEndB);
            
    if ((iStartA == iStartB) and (iEndA == iEndB) )
    {
        isTriangular=true;        

    }
    else
    {
        isTriangular=false;

        int lowerSet= iStartA < iStartB ? 0 : 1;

        if (lowerSet == 0)
        {
            assert(iEndA<iStartB);
            
        }
        else
        {
            assert(iEndB<iStartA);
        }
    }

};

template<int dims>
bool twoBodyPotential<dims>::containedInSetA(int i1,int i2) const
{
    return (i1 >= iStartA) and (i2 <=iEndA);
}

template<int dims>
bool twoBodyPotential<dims>::containedInSetB(int i1,int i2) const
{
    return (i1 >= iStartA) and (i2 <=iEndB);
};

template<int dims>
void twoBodyPotential<dims>::setPeriodic(std::array<Real,dims> lBox   )
{
    periodic=true;
    for(int d=0;d<dimensions;d++)
    {
        _lBox[d]=lBox[d];
    }

};





template class twoBodyPotential<3>;


