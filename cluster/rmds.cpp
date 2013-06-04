#include "rmds.h"
#include "tools.h"
#include "lbfgs.h"

Rmds::Rmds()
{
    this->_method = "RMDS";
    this->_basisRows = 1;
    this->_sizeOfNew = 1;
    this->_basisIndexMethod = 1;

    this->_parameters = new char[100]; // Reikia pirma inicializuoti kad galetume kopijuoti duomenis
    this->SetParameters();
}

/**
 * Apskaiciuoja santykinę projekciją
 *
 * Naudojasi keletu pasleptų parametrų
 * Realative taškų inicializavimas parenkama artimiausia bazinio vektoriaus projekcija
 * Quazi newtone naudojame max 200 iteracijų
 *
 * Galima projekcija yra tik 2
 * Algoritmas tik po viena taska dabar atideda
 *
 */
double* Rmds::Train()
{
    // Algoritmas tik po viena taska dabar atideda
    //this->_sizeOfNew = 7;
    int i, j, iTmp, iTmp2;

    if (this->_rows == 0 || this->_cols == 0) {
        this->LoadData(this->_pFileName);
    }

    if (this->_rows > 0 && this->_cols > 0) {
        int lbfgsParam; // parameter for lbfgs algortithm, depend on projection dimension
        if (this->_projDim < 3) {
            lbfgsParam = 2;
        } else {
            lbfgsParam = 3;
        }
        int* index = (int*) malloc (this->_rows * sizeof(int));

        double* pBasisData = (double*) malloc(this->_basisRows * this->_cols * sizeof(double));
        double* pBasisInitProj = (double*) malloc(this->_basisRows * this->_projDim * sizeof(double));
        double* pBasisProj = (double*) malloc(this->_basisRows * this->_projDim * sizeof(double));

        this->_pProj = (double*) realloc(this->_pProj, this->_rows * this->_projDim * sizeof(double));

        double* pNewData = (double*) malloc((this->_rows - this->_basisRows) * this->_cols * sizeof(double));

        this->getIndexes(index, this->_basisIndexMethod);

        // Baziniai vektoriai
        double* pY = pBasisData; // Pagreitinimui
        for (i = 0; i < this->_basisRows; i++) {
            iTmp = index[i] * this->_cols;
            for (j = 0; j < this->_cols; j++, pY++) {
                //*(pBasisData + i * this->_cols + j) = *(this->_pData + index[i] * this->_cols + j);
                *(pY) = *(this->_pData + iTmp + j);
            }
        }

        // Like vektoriai
        pY = pNewData; // Pagreitinimui
        for (i = this->_basisRows; i < this->_rows; i++) {
            iTmp = index[i] * this->_cols;
            for (j = 0; j < this->_cols; j++, pY++) {
                //*(pNewData + i * this->_cols + j) = *(this->_pData + index[i + this->_basisRows] * this->_cols + j);
                *(pY) =  *(this->_pData + iTmp + j);
            }
        }

        // rezultatai
        double* pResults = (double*) malloc(this->_rows * this->_projDim * sizeof(double));

        // Atvaizduojame bazinius vektorius i plokstuma
        this->InitProjectionVectors(pBasisData, this->_initMethod, this->_basisRows, this->_cols, this->_projDim, pBasisInitProj);
        this->Smacof(pBasisData, this->_basisRows, this->_cols, this->_iterations, this->_projDim, pBasisInitProj, pBasisProj);

        memcpy(pResults, pBasisProj, this->_basisRows * this->_projDim * sizeof(double));

        // Quazi-Newton
        ap::real_1d_array x; //, gr;
        //x.setbounds(1, this->_sizeOfNew * this->_projDim);
        double *pX, dMin;
        int info, k, rowsTaken;

        double* pBasisNewD = (double *) malloc(this->_basisRows * this->_sizeOfNew * sizeof(double));
        double* pNewNewD = (double *) malloc(this->_sizeOfNew * this->_sizeOfNew * sizeof(double));
        pX = (double *) malloc(this->_sizeOfNew * this->_cols * sizeof(double));

        pY = pResults + this->_basisRows * this->_projDim;

        for (int iRow = 0; iRow < (this->_rows - this->_basisRows); iRow++) {
            // Paruosiame reikiamas matricas relative gradiento skaiciavimui
            for (rowsTaken = 0; rowsTaken < this->_sizeOfNew && ((iRow + rowsTaken) < (this->_rows - this->_basisRows)); rowsTaken++) {
                iTmp = rowsTaken * this->_cols;
                iTmp2 = (iRow + rowsTaken) * this->_cols;
                for (j = 0; j < this->_cols; j++) {
                    *(pX + iTmp + j) = *(pNewData + iTmp2 + j);
                }
            }
            /*
            for (j = 0; j < this->_cols; j++) {
                *(pX + j) = *(pNewData + i * this->_cols + j);
            }*/
            Tools::EucDist2Matrix(pBasisData, pX,  this->_basisRows, rowsTaken, this->_cols, pBasisNewD);

            Tools::EucDist(pX, rowsTaken, this->_cols, pNewNewD);

            // Inicalizuojant pradiniu tasku priskiriame artimiausią bazinio vekrotoriaus projekciją
            x.setbounds(1, rowsTaken * this->_projDim);
            for (i = 0; i < rowsTaken; i++) {
                iTmp = i * this->_basisRows;
                dMin = *(pBasisNewD + iTmp); // pBasisNewD[0]
                k = i * this->_projDim;
                for (j = 1; j < this->_basisRows; j++) {
                    if (*(pBasisNewD + iTmp + j) < dMin) {
                        dMin = *(pBasisNewD + iTmp + j);
                        k = i * this->_projDim + j;
                    }
                }

                //x(1) = *(pBasisProj + k);
                //x(2) = *(pBasisProj + k + 1);
                for (j = 0; j < this->_projDim; j++) {
                    x(i * this->_projDim + j + 1) = *(pBasisProj + k + j);
                }
            }

            // Apskaiciuoti naujo tasko projekcija
            // (projDim, 3<=z<=7 or <=projDim, x, minGr, minF, minX, maxIter, info)
            lbfgsminimize(this->_projDim, lbfgsParam, x, 1E-8, 1E-8, 1E-8, 200, info, pBasisNewD, pNewNewD, pBasisProj, this->_basisRows, rowsTaken, this->_projDim);


            // Saugoti i masyva 2 mati
            for (j = 1; j <= rowsTaken * this->_projDim; j++) {
                *(pY++) = x(j);
            }

            if (rowsTaken > 1) {
                iRow += rowsTaken - 1;
            }

            //*(pY++) = x(1);
            //*(pY++) = x(2);
            //*(pResults + this->_projDim * (this->_basisRows + i)) = x(1);
            //*(pResults + this->_projDim * (this->_basisRows + i) + 1) = x(2);


        }

        // Atstatome projekciją
        pY = pResults;
        for (i = 0; i < this->_rows; i++) {
            iTmp = index[i] * this->_projDim;
            for (j = 0; j < this->_projDim; j++, pY++) {
                //*(this->_pProj + index[i] * this->_projDim + j) = *(pResults + this->_projDim * i + j);
                *(this->_pProj + iTmp + j) = *(pY);
            }
        }

        free(index);
        free(pBasisInitProj);
        free(pBasisProj);
        free(pBasisData);
        free(pNewData);
        free(pResults);

        return this->_pProj;
    }
}


/**
 * Funtcion set index array according selected type
 * 1 - random, 2 - pca with desc order
 *
 */
void Rmds::getIndexes(int* index, int type = 1)
{
    double* pInitProj;
    int i, rest, mod, t;

    switch (type) {
        //On line according PCA but
        // order is like this: 1, n, 2, n-1, 3, n-2, ...
        case 4:
            pInitProj = (double*) malloc(this->_rows * sizeof(double));
            this->InitProjectionVectors(this->_pData, 4, this->_rows, this->_cols, 1, pInitProj);
            Tools::BubbleSortDesc(pInitProj, this->_rows, index);
            int* index2;
            index2 = (int*) malloc(this->_rows * sizeof(int));
            t = 0;
            for (i = 0; i < this->_rows; i+=2) {
                index2[i] = index[t];
                if ((i + 1) < this->_rows) {
                    index2[i + 1] = index[this->_rows - 1 - t];
                }
                t++;
            }
            for (i = 0; i < this->_rows; i++) {
                *(index + i) = *(index2 + i);
            }
            free(index2);

        // On line according the biggest variance
        case 3:
            // Iskiriame atminti ir surandame didziausiu dispersiju projekcija duomenu i tiese
            pInitProj = (double*) malloc(this->_rows * sizeof(double));
            this->InitProjectionVectors(this->_pData, 3, this->_rows, this->_cols, 1, pInitProj);
            Tools::BubbleSortDesc(pInitProj, this->_rows, index);
            free(pInitProj);

            rest = this->_rows % this->_basisRows;
            mod = (int) (this->_rows - rest) / this->_basisRows;
            // Perkeliame vienodais tarpais isrinktus vektorius i pradzia
            // Likusius reiketu sumaisyti bet jei nepridesime ju prie baziniu
            // tai situacija nesikeicia
            for (i = 1; i < this->_basisRows; i++) {
                rest = index[i];
                index[i] = index[i * mod];
                index[i * mod] = rest;
            }
        // On line according PCA
        case 2:
            // Iskiriame atminti ir surandame pca projekcija duomenu i tiese
            pInitProj = (double*) malloc(this->_rows * sizeof(double));
            this->InitProjectionVectors(this->_pData, 4, this->_rows, this->_cols, 1, pInitProj);
            Tools::BubbleSortDesc(pInitProj, this->_rows, index);
            free(pInitProj);

            rest = this->_rows % this->_basisRows;
            mod = (int) (this->_rows - rest) / this->_basisRows;
            // Perkeliame vienodais tarpais isrinktus vektorius i pradzia
            // Likusius reiketu sumaisyti bet jei nepridesime ju prie baziniu
            // tai situacija nesikeicia
            for (i = 1; i < this->_basisRows; i++) {
                rest = index[i];
                index[i] = index[i * mod];
                index[i * mod] = rest;
            }
            break;
        // Random
        default:
        case 1:
            Tools::Shuffle(index, this->_rows, this->GetSeed());
            break;
    }
}
