#include <vector>

class cellData
{

};

class indexCellData : public cellData
{   
    public:

    using index_t = int;


    indexCellData():
    _nParticles(0)
    {}

    void resize(index_t newSize)  {_particleIndices.resize(newSize) ;
    
    }

    size_t size() const {return _nParticles;}


    index_t add(index_t iParticle){
        if (_nParticles >=_particleIndices.size())
        {
            _particleIndices.resize(_nParticles + buffer);
        }
        _particleIndices[_nParticles]=iParticle;
        _nParticles++;
        return _nParticles-1;
    }

    auto nParticles() const  {return _nParticles;}

    void remove(index_t subCellIndex)
    {
        std::swap(_particleIndices[_nParticles-1] , _particleIndices[subCellIndex]);
        _nParticles--;
    }

    void moveTo(int subCellIndex, indexCellData & cellDestination)
    {
        remove(subCellIndex);
        cellDestination.add( getParticleIndex(subCellIndex ));
    }  

    index_t  getParticleIndex(index_t i) const {
        return _particleIndices[i]; 
        } // get the i_th particle index

    private:  

    const int buffer = 10;
    std::vector<index_t> _particleIndices;
    index_t _nParticles; // number of particles containted in the box
};