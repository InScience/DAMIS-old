#include "dma.h"

Dma::Dma()
{
    this->_method = "DMA";
    this->_neighbour = 0;
    this->_parameters = new char[100];
    sprintf(this->_parameters, "k-factor:%d", _neighbour);
}

double* Dma::Train(bool zeidel = false)
{
    if (this->LoadData(this->_pFileName) > 0) {
        if (this->_rows>0 && this->_cols>0) {
            // projekcijos vektoriai
            this->_pProj = (double*) realloc(this->_pProj, this->_rows * this->_projDim * sizeof(double));
            double *pInitProj = (double*) malloc(this->_rows * this->_projDim * sizeof(double));

            srand((unsigned) this->GetSeed());
            if (this->InitProjectionVectors(this->_pData, this->_initMethod, this->_rows, this->_cols, this->_projDim, pInitProj) > 0) {
                if (zeidel) {
                    this->_method = "DMAZ";
                    this->EvaluateZ(this->_pData, this->_rows, this->_cols, this->_neighbour,this->_iterations, this->_projDim, pInitProj, this->_pProj);
                } else {
                    this->_method = "DMA";
                    this->Evaluate(this->_pData, this->_rows, this->_cols, this->_neighbour, this->_iterations, this->_projDim, pInitProj, this->_pProj);
                }
            }
           free(pInitProj);
           return this->_pProj;
        }
    }
    return NULL;
}

void Dma::Evaluate(double *pData, int rows, int cols, int K, int iterations, int projDim, double *initProj, double *proj)
{
    double *d, *GT_B, theSum, dist;
    int i, j, k, t, iter, pos;
    int *V;
    double *pGT_B, *pD, *pInitProj, *pProj, Y;

    /* Inicializuojame svoriu matrica */
    V = (int*) calloc((rows), sizeof(int));

    /* create pointers to distance and GT_B matrices */
    /* Setting matrices with values of zero */
    d = (double*) calloc((rows * rows), sizeof(double));
    GT_B = (double*) calloc((rows * rows), sizeof(double));

    /* Matrices of distances */
    Tools::EucDist(pData, rows, cols, K, d);

    //memcpy(
    for (i = 0; i < rows * projDim; i++) {
        *(proj + i) = *(initProj + i);
    }

    /* Iterational process */
    for (iter = 1; iter <= iterations; iter++) {
        /* Creating Guttman tranformation matrix GT_B*/
        /* We changed original GT_B with GT_B = (GT_B - V) */
        pos = -rows;
        for (i = 0; i < rows; i++) {
            k = i - K;
            if (i - K < 0) {
                k = 0;
            }
            theSum = 0.0;
            pos += rows;
            pGT_B = GT_B + pos + k; // Pointer into GT_B element (row)

            /* Elements before diagonal */
            for (j = k; j < i; j++, pGT_B++) {
                theSum -= *(pGT_B);
            }

            pGT_B++;

            pD = d + pos + i + 1;

            /* Elements after diagonal */
            for (j = i + 1; (j < rows && j <= i + K); j++, pGT_B++, pD++) {
                dist = 0.0;
                for (t = 0; t < projDim; t++) {
                    // Best performance versus pointes or dist = Y * Y - vistik > 4% pagreitejimas
                    Y = *(initProj + projDim * i + t) - *(initProj + projDim * j + t);
                    dist += Y * Y;
                    //dist += (*(initProj + projDim * i + t) - *(initProj + projDim * j + t)) * (*(initProj + projDim * i + t) - *(initProj + projDim * j + t));
                }
                dist = sqrt(dist);

                if (dist > 0.000001) {
                     *(pGT_B) = 1.0 - *(pD) / dist;   // For saving space of V matrix
                } else {                             // we are adding V elements direct to GT_B
                    *(pGT_B) = 1.0;
                }
//                 *(pGT_B) += 1;
                *(GT_B + rows * j + i) = *(pGT_B);
                theSum -= *(pGT_B);
            }
            *(GT_B + pos+ i) = theSum;
            *(V + i) = j - k - 1;  // Diagonal elements of matrix V
        }

        /* Calculating new projection of points */
        // proj = proj + 0.5 * diag(diag(V).^(-1)) * (GT_B - V) * proj;
        pos = -rows;
        pProj = proj;
        pInitProj = initProj;
        for (i = 0; i < rows; i++) {
            k = i - K;
            if (i - K < 0) {
                k = 0;
            }

            pos = pos + rows;
            /* Projection coords */
            for (t = 0; t < projDim; t++, pProj++) {
                theSum = 0.0;
                pInitProj = initProj + projDim * k + t;
                for (j = k; (j < rows && j <= i + K); j++) {
                    theSum += *(GT_B + pos + j) * *(pInitProj);
                    pInitProj += projDim;
                }
                *(pProj) += theSum * 0.5 / *(V + i);
            }
        }

        for (i = 0; i < rows * projDim; i++) {
            *(initProj + i) = *(proj + i);
        }
    }
}

/************************************************************/
/* After each iteration this algoritm change point order
 * this whay lets us to receive better error */
void Dma::EvaluateZ(double *pData, int rows, int cols, int K, int iterations, int projDim, double *initProj, double *proj)
{
    double *d, *GT_B, theSum, dist;
    int i, j, k, t, iter, pos, *index, rowI, rowJ;
    int *V;
    double *pGT_B, Y;

    /* Inicializuojame svoriu matrica */
    V = (int*) calloc((rows), sizeof(int));

    /* create pointers to distance and GT_B matrices */
    /* Setting matrices with values of zero */
    d = (double*) calloc((rows * rows), sizeof(double));
    GT_B = (double*) calloc((rows * rows), sizeof(double));
    index = (int*) calloc(rows, sizeof(int));


    /* Matrices of distances */
    Tools::EucDist(pData, rows, cols, d);

    for (i = 0; i < rows * projDim; i++) {
        *(proj + i) = *(initProj + i);
    }

    /* Iterational process */
    for (iter = 1; iter <= iterations; iter++) {
        /* Creating Guttman tranformation matrix GT_B*/
        /* We changed original GT_B with GT_B = (GT_B - V) */
        Tools::Shuffle(index, rows, this->GetSeed());
        pos = -rows;
        for (i = 0; i < rows; i++) {
            rowI = index[i];
            k = i - K;
            if (i - K < 0)
                k = 0;
            theSum = 0.0;
            pos += rows;
            pGT_B = GT_B + pos + k; // Pointer into GT_B element (row)

            /* Elements before diagonal */
            for (j = k; j < i; j++, pGT_B++) {
                theSum -= *(pGT_B);
            }

            pGT_B++;

            /* Elements after diagonal */
            for (j = i + 1; (j < rows && j <= i + K); j++, pGT_B++) {
                rowJ = index[j];
                dist = 0.0;
                for (t = 0; t < projDim; t++) {
                        // Best performance versus pointes or dist = Y * Y
                        //Y = *(proj + projDim * rowI + t) - *(proj + projDim * rowJ + t);
                        //dist += Y * Y;
                        // Kad eksperimento salygos butu vienodos
                        dist += (*(proj + projDim * rowI + t) - *(proj + projDim * rowJ + t)) * (*(proj + projDim * rowI + t) - *(proj + projDim * rowJ + t));
                }
                dist = sqrt(dist);

                if (dist > 0.000001) {
                     *(pGT_B) = 1.0 - *(d + rows * rowI + rowJ) / dist;   // For saving space of V matrix
                } else {                             // we are adding V elements direct to GT_B
                    *(pGT_B) = 1.0;
                }
                 *(GT_B + rows * j + i) = *(pGT_B);
                theSum -= *(pGT_B);
            }
            *(GT_B + pos + i) = theSum;
            *(V + i) = j - k - 1;  // Diagonal elements of matrix V
        }

        /* Calculating new projection of points */
        // proj = proj + 0.5 * diag(diag(V).^(-1)) * (GT_B - V) * proj;
        pos = -rows;
        for (i = 0; i < rows; i++) {
            rowI = index[i];
            k = i - K;
            if (i - K < 0)
                k = 0;

            pos = pos + rows;
            /* Projection coords */
            for (t = 0; t < projDim; t++) {
                theSum = 0.0;
                for (j = k; (j < rows && j <= i + K); j++) {
                    theSum += *(GT_B + pos + j) * *(initProj + projDim * index[j] + t);
                }
                *(proj + projDim * rowI + t) += theSum * 0.5 / *(V + i);
            }
        }

        for (i = 0; i < rows * projDim; i++) {
            *(initProj + i) = *(proj + i);
        }
    }
}
