
using Real = double;

template<int dims,class V_t>
void addTwoBodyIsotropicForcesRectangular(
    const V_t & V,
    const Real * particles,
    Real * forces,
     int i1,int i2, // i-particle range
     int j1, int j2, // j-particle range
     int t0 , int t1, // time range
     int N, // total number of particles 
     int T, // total number of time slices
     const Real * lBox
     )
{
    /*
    Assumes no intersection between i-particle range and j-particle range 
    */
    Real sum2b=0;
    constexpr int D = dims;

    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for schedule(static) collapse(3)
    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<=j2;jParticle++)
            {
                Real r2=0;
                std::array<Real,3> diff;

                for(int d=0;d<dims;d++)
                {
                    diff[d]=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diff[d]= utils::differencePBC(diff[d],lBox[d],lBoxInverse[d]);
                    r2+= diff[d] * diff[d];
                }

                auto r = std::sqrt(r2);
                auto rInverse = 1./r;
                auto dVdr = V.radialDerivative(r);                
                
                for(int d=0;d<dims;d++)
                {
                    forces[iParticle + d*N + t*N*D ]+=dVdr*diff[d]*rInverse;
                    forces[jParticle + d*N + t*N*D ]-=dVdr*diff[d]*rInverse;
                }
             }
        }
}



template<int dims,class V_t>
void addTwoBodyIsotropicForcesTriangular(
    const V_t & V,
    const Real * particles,
    Real * forces,
     int i1,int i2, // i-particle range
     int j1, // beginning of the i-range
     int t0 , int t1, // time range
     int N, // total number of particles 
     int T, // total number of time slices
     const Real * lBox
     )
{
    /*
    Assumes no intersection between i-particle range and j-particle range 
    */
    Real sum2b=0;
    constexpr int D = dims;

    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }


    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<iParticle;jParticle++)
            {
                Real r2=0;
                std::array<Real,3> diff;

                for(int d=0;d<dims;d++)
                {
                    diff[d]=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diff[d]= utils::differencePBC(diff[d],lBox[d],lBoxInverse[d]);
                    r2+= diff[d] * diff[d];
                }

                auto r = std::sqrt(r2);
                auto rInverse = 1./r;
                auto dVdr = V.radialDerivative(r);                
                
                for(int d=0;d<dims;d++)
                {
                    forces[iParticle + d*N + t*N*D ]+=dVdr*diff[d]*rInverse;
                    forces[jParticle + d*N + t*N*D ]-=dVdr*diff[d]*rInverse;
                }
             }
        }
}







template<int dims,class V_t>
Real evaluateTwoBodyRectangular(
    const V_t & V,
    const Real * particles,
     int i1,int i2, // i-particle range
     int j1, int j2, // j-particle range
     int t0 , int t1, // time range
     int N, // total number of particles 
     int T, // total number of time slices
     const Real * lBox
     )
{
    Real sum2b=0;
    constexpr int D = dims;


    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(3)
    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<=j2;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<dims;d++)
                {
                    Real diffd=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diffd= utils::differencePBC(diffd,lBox[d],lBoxInverse[d]);
                    r2+= diffd * diffd;
                }

                sum2b+=V(std::sqrt(r2));
        }
    }

return sum2b;
}

template<int dims,class V_t>
Real evaluateTwoBodyTriangular(
    const V_t & V,
    const Real * particles,
     int i1,int i2, // i-particle range
     int j1 ,// beginning of the i-particle set
     int t0 , int t1, // time range
     int N, // total number of particles 
     int T, // total number of time slices
     const Real * lBox
     )
{
    Real sum2b=0;
    constexpr int D = dims;


    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(2)
    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<iParticle;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<dims;d++)
                {
                    Real diffd=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diffd= utils::differencePBC(diffd,lBox[d],lBoxInverse[d]);
                    r2+= diffd * diffd;
                }

                sum2b+=V(std::sqrt(r2));
        }
    }

return sum2b;
}



template<int dims>
template<class V_t>
Real twoBodyPotential<dims>::operator()(const V_t & V, // two body potential
const Real * positions, // raw data in a contigous array of shape (N,dims, T)
int i1 , int i2, // particle range updated 
 int t0 , int t1,  // time range updated
 int N, int T // shape information of input data
 ) const
{

    if ( not isTriangular )
    {
        bool isInA = containedInSetA(i1,i2);
        bool isInB = containedInSetB(i1,i2);
        
        if (isInA)
        {
            return evaluateTwoBodyRectangular<dims,V_t>(V,positions,i1,i2,iStartB,iEndB,t0 , t1,N,T,_lBox.data());
        }
        else if (isInB)
        {
           return evaluateTwoBodyRectangular<dims,V_t>(V,positions,i1,i2,iStartA,iEndA,t0 , t1,N,T,_lBox.data());
        }

        return 0;
    }
    else
    {
        bool isInA = containedInSetA(i1,i2);


        if (isInA)
        {
            auto sum = evaluateTwoBodyTriangular<dims,V_t>(V,positions,i1,i2,iStartA,t0 , t1,N,T,_lBox.data());
            sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,i2+1,iEndA,i1,i2,t0 , t1,N,T,_lBox.data() );


            return sum;
        }
        else
        {
            return 0;
        }
        
    }
}


template<int dims>
template<class V_t>
void twoBodyPotential<dims>::addForce(
const V_t & V, // two body potential
const Real * positions, // raw data in a contigous array of shape (N,dims, T),
 Real * forces, // raw data in a contigous array of shape (N,dims, T) where to add the force for each particle
int i1 , int i2, // particle range updated 
 int t0 , int t1,  // time range updated
 int N, int T // shape information of input data
 ) const
{
    if ( not isTriangular )
    {
        bool isInA = containedInSetA(i1,i2);
        bool isInB = containedInSetB(i1,i2);
        
        if (isInA)
        {
            addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,i1,i2,iStartB,iEndB,t0 , t1,N,T,_lBox.data());
        }
        else if (isInB)
        {
           addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,i1,i2,iStartB,iEndB,t0 , t1,N,T,_lBox.data());
        }

    }
    else
    {
        bool isInA = containedInSetA(i1,i2);

        if (isInA)
        {
            addTwoBodyIsotropicForcesTriangular<dims,V_t>(V,positions,forces,i1,i2,iStartA,t0 , t1,N,T,_lBox.data());
             addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,i2+1,iEndA,i1,i2,t0 , t1,N,T,_lBox.data() );
        }
        else
        {
            
        }
        
    }
}