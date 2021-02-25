#include "particleData.h"
#include <iostream>
#include "lattice.h"
#include <random>
#include "utils.h"
#include "twoBodyDistances.hpp"
#include "timers.h"



template<class testParticleData_t>
void timeTwoBodyPotential(std::string method, size_t N, int nTrials)
{

    ADD_TIMER("2bDistances");


     std::cout << "Timing particle data containers" << std::endl;
    
     Real cut_off=1;
    constexpr int D = 3;
    Real density = 1e-3;

    Real l=std::pow(N/density,1./3);


    std::array<Real, D> lBox={l,l,l};
    std::array<size_t,D> nBoxes;
    std::array<Real,D> leftEdge={-lBox[0]/2.,-lBox[1]/2.,-lBox[2]/2.};
    std::array<Real,D> rightEdge={lBox[0]/2.,lBox[1]/2.,lBox[2]/2.};

    for (int d=0;d<D;d++)
    {
        nBoxes[d]=std::min(std::pow(N,1./3),std::floor(lBox[d]/(cut_off) ) );
    }



    auto lattptr= std::shared_ptr<lattice>( new  lattice(nBoxes, leftEdge,rightEdge));


    //ASSERT_EQ(lattptr->lowIndex(0),1);
    //ASSERT_EQ(lattptr->highIndex(0),10);

    lattptr->checkNeighbourIndexing();


    double * positions = (Real *)aligned_alloc( 128, N*D*sizeof(Real)  );
    double * positionsPBC = (Real *)aligned_alloc( 128, N*D*sizeof(Real)  );


    std::default_random_engine generator(12);


    utils::initRandom(positions,3,0,N-1,N,generator,-l/2,l/2);


    utils::restrictToBox(positionsPBC,positions,0,N-1,N,leftEdge,lBox);


    utils::gaussianInteraction gauss(2,1);
    utils::tetaInteraction teta(1,1);



    if (method == "direct")
    {
        Real vCheck = utils::twoBodyDistancesIsotropicReduction<utils::gaussianInteraction,3>(positions,gauss,N,lBox);
        std::cout << "Sum: " << vCheck << std::endl;

        for (int n=0;n<nTrials;n++)
    {

        int iStart=0;
        int iEnd=N-1;

       
        utils::initRandom(positions,3,iStart,iEnd,N,generator,-l/2,l/2);


        START_TIMER("2bDistances") ;
        Real vCheckUpdate = utils::twoBodyDistancesIsotropicReduction<utils::gaussianInteraction,3>(positions,gauss,N,lBox);
        STOP_TIMER("2bDistances") ;


        std::cout << vCheckUpdate << std::endl;
        

    }

    }


    else{

        testParticleData_t particleAcc(lattptr);
        particleAcc.list(positionsPBC,0,N-1,N);
        auto v = twoBodyDistancesIsotropicReduction(particleAcc,gauss,positionsPBC,0,N-1,N);
    
    
    for (int n=0;n<nTrials;n++)
    {

        int iStart=0;
        int iEnd=N-1;

        utils::initRandom(positions,3,iStart,iEnd,N,generator,-l/2,l/2);
         START_TIMER("2bDistances") ;
        utils::restrictToBox(positionsPBC,positions,iStart,iEnd,N,leftEdge,lBox);
        particleAcc.updateList(positionsPBC,iStart,iEnd,N);
        Real vUpdate=twoBodyDistancesIsotropicReduction(particleAcc,gauss,positionsPBC,iStart,iEnd,N);
        STOP_TIMER("2bDistances") ;
        
        std::cout << vUpdate << std::endl;
    }

    }

    free(positions);
    free(positionsPBC);

    std::cout << timers::getInstance().report() << std::endl;

}


int main(int argc,char** argv)
{
    int nSamples=10;
    int N=100;
    std::string method="direct";
    
    if ( argc > 1) method = argv[1];
    if ( argc > 2) N = atoi( argv[2] );
    if (argc > 3 ) nSamples=atoi( argv[3] );

    

    if (method == "direct" or "acc_index")
        {
            timeTwoBodyPotential<particleDataIndex>(method,N,nSamples);
        }
        else if (method == "acc_copy")
        {
            timeTwoBodyPotential<particleData3D>(method,N,nSamples);
        }
}