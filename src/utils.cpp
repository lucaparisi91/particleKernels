#include "utils.h"
#include <cmath>
#include <cassert>
namespace utils
{
    Real differencePBC(Real t, Real lBox, Real lBoxInverse ) {return ( t - std::round(t*lBoxInverse )*lBox);}

    Real restrictToBox(Real x, Real left, Real lBox, Real lBoxInverse)
    {

        x= x - std::round( ((x-left)*lBoxInverse - 0.5 ))*lBox;
        assert(x>=left);
        assert(x<=left+lBox);
        
        return x;
    };

    void restrictToBox(Real * positionsPBC, const Real * positionsOld, int iStart,int iEnd, int N, const std::array<Real,3> & left, const std::array<Real,3> & lBox )
    {
        std::array<Real,3> lBoxInverse{1./lBox[0],1./lBox[1],1./lBox[2]};
        for(int i=iStart;i<=iEnd;i++)
        {
            for(int d=0;d<3;d++)
            {
                positionsPBC[i+d*N]=restrictToBox(positionsOld[i+d*N],left[d],lBox[d],lBoxInverse[d]);
            }
        }
    }


};
