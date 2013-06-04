#ifndef MDS_H
#define MDS_H

#include <math.h> // sqrt
#include <iostream>
#include <cassert>
#include <iomanip> // setw, setprecision
#include "tools.h"
#include <string.h> // arba cstring
#include "libs/alglib/ap.h"
#include "alglib_pca.h"

using namespace std;

class Mds
{
    public:
        Mds();
        ~Mds();

        void SetFile(char* pFileName) { _pFileName = pFileName; }
        void SetIterations(char* pIterations) { _iterations = atoi(pIterations); }
        void SetInitMethod(char* pInitMethod) { _initMethod = atoi(pInitMethod); }
        void SetRank(int rank) { _rank = rank; }
        void SetSeed(int seed) { _seed = seed; }
        int GetSeed() { return _seed; }
        double* GetProjection() { return _pProj; }
        char* GetMethod() { return _method; };
        char* GetParameters() { return _parameters; };

        double* Train(bool zeidel = false);
        double Error(int type);
        double MinimalWiringError();
        double BadClassifiedVectorsError();
        double SpearmanCoefficient();
        void SaveResults(char* pFile, double time);
        void SaveResults(char* pFile, double time, double error, double* pData);

        int LoadData(char *pFileName);
        int InitProjectionVectors(double* pData, int method, int rows, int cols, int projDim, double* pPlaneArr);
        void Smacof(double *pData, int rows, int cols, int iterations, int projDim, double *initProj, double *proj);
        void SmacofZ(double *pData, int rows, int cols, int iterations, int projDim, double *initProj, double *proj);

        char* _method;
        char* _parameters;
        int _rank;
        char* _pFileName;
        int _iterations;
        int _initMethod;
        int _projDim;
        int _rows;
        int _cols;
        int _seed;
        double* _pData;
        double* _pProj;
    protected:

        
    private:
        

};

#endif
