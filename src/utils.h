using Real = double;
#include <array>


namespace utils
{
    Real differencePBC(Real t, Real lBox, Real lBoxInverse ) ;  
    Real restrictToBox(Real x, Real left, Real lBox, Real lBoxInverse);

    void restrictToBox(Real * positionsPBC, const Real * positionsOld, int iStart,int iEnd, int N, const std::array<Real,3> & left, const std::array<Real,3> & lBox );
    

};