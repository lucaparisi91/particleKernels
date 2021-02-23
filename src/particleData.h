#include "lattice.h"
#include "cellData.h"
#include "cassert"

template<class T>
class particleDataBase
{
    public:

    using cell_t = T;
    using index_t = int;

    particleDataBase(std::shared_ptr<lattice> _lattice);
    
    virtual void list(const Real * positions, int iStart, int iEnd, int N)=0;



    //virtual void updateList(const Real * positions, int iStart,int iEnd, int N)=0;

    const auto & getLattice() const {return *(_lattice);}


    int cellIndex(int iParticle) const  {return _cellIndexPerParticle[iParticle];}

    int subCellIndex(int iParticle) const {
        
        return _subCellIndexPerParticle[iParticle];        
     }

    auto & operator[](size_t i)  {return *_cellData[i];}
    const auto & operator[](size_t i) const {return *_cellData[i];}


    void recordParticle(index_t iParticle,index_t iCell, index_t iSubCell) {
        _cellIndexPerParticle[iParticle]=iCell;
        _cellIndexPerParticle[iSubCell]=iSubCell;

    }

    void resizeParticleData(size_t nMax)
    {
        _cellIndexPerParticle.resize(nMax,-1);
        _subCellIndexPerParticle.resize(nMax,-1);
    }


    void recordParticleCellPosition(index_t iParticle,index_t iCell , index_t iSubCell)
    {
        _cellIndexPerParticle[iParticle]=iCell;
        _subCellIndexPerParticle[iParticle]=iSubCell;
    }

    void checkPBCGhostCells();

    private:


    std::shared_ptr<lattice> _lattice;

    /* Per particle data */
    std::vector<int>  _cellIndexPerParticle;
    std::vector<int> _subCellIndexPerParticle;

    /* Cell data  */
    std::vector<std::shared_ptr<cell_t> > _cellData;
};


class particleDataIndex : public particleDataBase<indexCellData>
{
    public:

    particleDataIndex(std::shared_ptr<lattice> lattice) : particleDataBase<indexCellData>(lattice) {}



    void list(const double * positions, int iStart, int iEnd, int N);

    private:

};


namespace utils
{
    Real differencePBC(Real t, Real lBox, Real lBoxInverse );

};
template<class V_t>
Real twoBodyDistancesIsotropicReduction( const particleDataIndex & container , const V_t & op, const Real * particles , int iStart, int iEnd , int N   )
{
    Real sum2b=0;

    


    const auto & currentLattice = container.getLattice();
    const auto dimensions = currentLattice.dimensions();

    // loop over particles to be updated
    for (int iParticle=iStart;iParticle<=iEnd;iParticle++)
    {
        int iCell=container.cellIndex(iParticle);
        int iSubCell=container.subCellIndex(iParticle);

        assert(iCell<currentLattice.extendedSize());
        assert(iSubCell<container[iCell].size());


        for (int j=0;j<currentLattice.nCellsNeighbourhood() ;j++)
        {
            int jCell=currentLattice.getNeighbour(iCell,j);
            assert(jCell<currentLattice.extendedSize());
            

            for (int jj=0;jj<container[jCell].nParticles();jj++)
            {
                int jParticle=container[jCell].getParticleIndex(jj);

                assert(jParticle<N);

                if ( jParticle < iParticle)
                {
                    Real r2=0;
                    for (int d=0;d<dimensions;d++)
                    {

                        Real diffd=( particles[iParticle + d*N ] - (particles[jParticle + d*N] 
                        + currentLattice.wrap(jCell,d) 
                        )
                           );

                        //diffd=utils::differencePBC(diffd,currentLattice.lengthBox()[d],1./currentLattice.lengthBox()[d] );
                        

                        
                        r2+=diffd*diffd;
                    }
                    sum2b+=op(std::sqrt(r2));
                }

            }

        }

    }
    return sum2b;
}