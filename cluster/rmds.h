#ifndef RMDS_H
#define RMSD_H

#include "mds.h"
#include "tools.h"
#include "lbfgs.h"

using namespace std;

class Rmds : public Mds
{
    public:
        Rmds();
        //~Rmds();

        void SetBasisRows(char* pBasisRows)
        {
            _basisRows = atoi(pBasisRows);
            this->SetParameters();
	};
        void SetSizeOfNew(char* pSizeOfNew) { _sizeOfNew = atoi(pSizeOfNew); };
        void SetBasisIndextMethod(char* pBasisIndexMethod)
        {
            _basisIndexMethod = atoi(pBasisIndexMethod);
            this->SetParameters();
        };
        void SetParameters()
        {
            sprintf(this->_parameters, "basis:%d basisIndexMethod:%d sizeOfNew:%d", _basisRows, _basisIndexMethod, _sizeOfNew);
        }
        void getIndexes(int* index, int type);
        double Error(int type = 0) { return Mds::Error(type); };

        double* Train();

    protected:

    private:
        int _basisRows;
        int _sizeOfNew;
        int _basisIndexMethod;
};

#endif
