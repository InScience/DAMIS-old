#ifndef DMA_H
#define DMA_H 

#include "mds.h"

using namespace std;

class Dma : public Mds
{
    public:
        Dma();
        //~Dma();

        double* Train(bool zeidel);
        
        void SetNeighbour(char* pNeighbour)
        {
            _neighbour = atoi(pNeighbour);
            sprintf(this->_parameters, "k-factor:%d", _neighbour);
        };
        void Evaluate(double *pData, int rows, int cols, int K, int iterations, int projDim, double *initProj, double *proj);
        void EvaluateZ(double *pData, int rows, int cols, int K, int iterations, int projDim, double *initProj, double *proj);

    protected:


    private:
        int _neighbour;
};

#endif


