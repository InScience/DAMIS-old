#include "mds.h"

Mds::Mds()
{
    this->_method = "MDS";
    this->_parameters = "";
    this->_pFileName = NULL;
    this->_iterations = 0;
    this->_initMethod = 0;
    this->_projDim = 2;
    this->_rows = 0;
    this->_cols = 0;
    this->_seed = (int) time(NULL);
    this->_pData = (double*) malloc(1000 * 4 * sizeof(double));
    this->_pProj = (double*) malloc(1000 * 4 * sizeof(double));
}

Mds::~Mds()
{
    free(this->_pData);
    free(this->_pProj);
}

/**
 *
 *
 * bool zeidel - if it is true than run zeidel type modification
 */
double* Mds::Train(bool zeidel)
{
    if (this->LoadData(this->_pFileName) > 0) {
        if (this->_rows>0 && this->_cols>0) {
            this->_pProj = (double*) realloc(this->_pProj, this->_rows * this->_projDim * sizeof(double));
            double *pInitProj = (double*) malloc(this->_rows * this->_projDim * sizeof(double));

            if (this->InitProjectionVectors(this->_pData, this->_initMethod, this->_rows, this->_cols, this->_projDim, pInitProj) > 0) {
                if (zeidel) {
                    this->SmacofZ(this->_pData, this->_rows, this->_cols, this->_iterations, this->_projDim, pInitProj, this->_pProj);
                } else {
                    this->Smacof(this->_pData, this->_rows, this->_cols, this->_iterations, this->_projDim, pInitProj, this->_pProj);
                }
            }
            free(pInitProj);
            return this->_pProj;
        }
    }
    return NULL;
}
/**
 * Type 0 - MDS error
 */
double Mds::Error(int type = 0)
{
    double* pDist = (double*) malloc(this->_rows * this->_rows * sizeof(double));
    double *pProjDist = (double*) malloc(this->_rows * this->_rows * sizeof(double));

    Tools::EucDist(this->_pData, this->_rows, this->_cols, pDist);
    Tools::EucDist(this->_pProj, this->_rows, this->_projDim, pProjDist);

    double error = Tools::MdsError(pDist, pProjDist, this->_rows);
    free(pDist);
    free(pProjDist);

    return error;
}

double Mds::MinimalWiringError()
{
    return Tools::MinimalWiringError(_pData, _rows, _cols, _pProj, _projDim);
}

double Mds::BadClassifiedVectorsError()
{
    return Tools::BadClassifiedVectors(_pData, _rows, _cols, _pProj, _projDim);
}

double Mds::SpearmanCoefficient()
{
    return Tools::ApSpearmanCoefficient(_pData, _rows, _cols, _pProj, _projDim);
}

void Mds::SaveResults(char* pFile, double time)
{
    ofstream resultsFile;
    resultsFile.open (pFile);
    if (resultsFile.is_open()) {
        resultsFile << "[Method]:" << this->GetMethod() << "\n";
        resultsFile << "[InitMethod]:" << this->_initMethod << "\n";
        resultsFile << "[Parameters]:" << this->GetParameters() << "\n";
        resultsFile << "[DataSet]:" << this->_rows << "x" << this->_cols << "\n";
        resultsFile << "[T]:" << fixed << setprecision(3) << time << "\n";
        resultsFile << "[I]:" << this->_iterations << "\n";
        resultsFile << "[E]:" << fixed << setprecision(15) << this->Error() << "\n";
        resultsFile << "[MWE]:" << fixed << setprecision(15) << this->MinimalWiringError() << "\n";
        //resultsFile << "[BCVE]:" << fixed << setprecision(15) << this->BadClassifiedVectorsError() << "\n";
        resultsFile << "[SC]:" << fixed << setprecision(15) << this->SpearmanCoefficient() << "\n";
		//resultsFile << "[SC]:" << fixed << setprecision(15) << 0.0 << "\n";
		
        for (int i=0; i<this->_rows; i++) {
            for (int j=0; j<this->_projDim; j++) {
                resultsFile << *(_pProj + i * this->_projDim + j) << ";";
            }
            resultsFile <<"3;"<< endl;
        }
        resultsFile.close();
    } else {
        printf("negali atidaryti failo!\n");
    }
}

void Mds::SaveResults(char* pFile, double time, double error, double* pData)
{
    ofstream resultsFile;
    resultsFile.open (pFile);
    if (resultsFile.is_open()) {
        resultsFile << "[T]:" << fixed << setprecision(3) << time << "\n";
        resultsFile << "[I]:" << this->_iterations << "\n";
        resultsFile << "[E]:" << fixed << setprecision(15) << error << "\n";
        resultsFile << "[MWE]:" << fixed << setprecision(15) << this->MinimalWiringError() << "\n";
        resultsFile << "[BCVE]:" << fixed << setprecision(15) << this->BadClassifiedVectorsError() << "\n";
        resultsFile << "[SC]:" << fixed << setprecision(15) << this->SpearmanCoefficient() << "\n";

        for (int i=0; i<this->_rows; i++) {
            for (int j=0; j<this->_projDim; j++) {
                resultsFile << *(pData + i * this->_projDim + j) << ";";
            }
            resultsFile <<"3;"<< endl;
        }
        resultsFile.close();
    } else {
        printf("negali atidaryti failo!\n");
    }
}

// Nuskaito duomenis, suranda stulpeliu ir eiluciu skaiciu
// pFileName - kelias iki duomenu failo
int Mds::LoadData(char *pFileName)
{
    double num;
    int i, j;

    ifstream file;
    file.open(pFileName, ifstream::in);
    if (file.is_open()) {
        //Nunuliname pradinius masyvo matmenis
        _rows = -1;
        _cols = 1000;
        // Pasiimame masyvo dydi
        j = 0;
        while (!file.eof()) {
            i = 0;
            string line;
            getline(file, line, '\n');
            char *tmpstr, *tmpstr2;
            tmpstr = &line[0];
            while (strlen(tmpstr) > 0 && i<_cols) {
                num = double (strtod(tmpstr, &tmpstr2));
                if (tmpstr != tmpstr2) { // do new number is found
                    tmpstr = tmpstr2;
                    i++;
                } else {
                    break;
                }
            }

            // Pirma gerai nuskaityta eilute jei sulpeliu maziau negu 1 arba nesutampas su cols skaiciumi
            if (j == 0) {
                if (i > 1) { _cols = i; } else { j = -1; }
            } else {
                if (i != _cols) { j--; }
            }
            j++;
        }
        _rows = j;

        file.close();
    }

    if (_rows > 0 && _cols > 0) {
        // pakartotinai nuskaitome faila dabar jau zinodami masyvo dydi
        // isskiriam vietos duomenims
        _pData = (double*) realloc (_pData, _rows * _cols * sizeof(double));

        file.open(pFileName, ifstream::in);
        if (file.is_open()) {
            j = 0;
            while (!file.eof() && j < _rows) {
                i = 0;
                string line;
                getline(file, line, '\n');
                char *tmpstr, *tmpstr2;
                tmpstr = &line[0];
                while (strlen(tmpstr) > 0 && i < _cols) {
                    num = double(strtod(tmpstr, &tmpstr2));
                    *(_pData + j * _cols + i) = num;
                    tmpstr = tmpstr2;
                    i++;
                }

                // Pirma gerai nuskaityta eilute jei sulpeliu maziau negu 1 arba nesutampas su cols skaiciumi
                if (j == 0) {
                    if (i < _cols) {
                      j = -1;
                    }
                } else {
                    if (i != _cols) {
                        j--;
                    }
                }
                j++;
            }
            file.close();
            return 1;
        }
    }
    return 0;
}

/**
 * Priskiria pradiniu vektoriu reiksmes
 *
 * @param pData     data
 * @param method    1,2,3,4
 * @param rows      data rows
 * @param cols      data cols
 * @param projDim   projection dimension
 * @param pPlaneArr return inicalized projection array
 * @return int      0|1. 0 - error appear
 */
int Mds::InitProjectionVectors(double *pData, int method, int rows, int cols, int projDim, double *pPlaneArr)
{
    int i, k;
    switch (method) {
        default:
        // Uzpildymas ant tieses (1000x50 suskaičiuoja per 0,02s)
        case 1:
            for (i=0; i<rows; i++) {
                for (k=0; k<projDim; k++) {
                    *(pPlaneArr + i * projDim + k) = (double) (i + 1) * (0.1 - 1.0 * k) + (double)(10 * k + 1) / 3 + 0.5;
                }
            }
            return 1;
            break;

        // Atsitiktinis uzpildymas (1000x50 suskaičiuoja per 0,08s)
        case 2:
            srand((unsigned) this->GetSeed());
            for (i=0; i<rows; i++) {
                for (k=0; k<projDim; k++) {
                    *(pPlaneArr + i * projDim + k) = (double) rand() *  10 / RAND_MAX - 5;
                }
            }
            return 1;
               break;

        // Didziausios dispersijos (1000x50 suskaičiuoja per 0,70s)
        case 3:
            double *pVar, *X;
            int *index;

            pVar = (double *) calloc (cols, sizeof(double)); // iskirame dispersijoms vieta
            index = (int *) calloc (cols, sizeof(int));  // iskiriame vietą indeksams
            Tools::CalcColsVariance(pData, rows, cols, pVar); // skaiciuojam stulpeliu dispersijas

            Tools::BubbleSortDesc(pVar, cols, index);

            // Priskirame koordinates
            X = pPlaneArr;
            for (i=0; i<rows; i++) {
                for (k=0; k<projDim; k++, X++) {
                    //*(pPlaneArr + i * projDim + k) = *(pData + i * this->_cols + index[k]);
                    *(X) = *(pData + i * cols + index[k]);
                }
            }
            return 1;

            free(pVar);
            free(index);
            break;
        // PCA (1000x50 suskaičiuoja per 50,3s)
        case 4:
            // Performatuojame i alglib formata
            ap::real_2d_array apData;
            apData.setbounds(0, rows-1, 0, cols-1);
            for (i=0; i<rows; i++) {
                for (k=0; k<cols; k++) {
                    apData(i, k) = *(pData + i * cols + k);
                }
            }
            calculatePcaProjection(apData, rows, cols, projDim, pPlaneArr);

            return 1;
            break;
    }
    return 0;
}

void Mds::Smacof(double *pData, int rows, int cols, int iterations, int projDim, double *initProj, double *proj)
{
    double *d, *GT_B, theSum, dist, Y;
    int i, j, t, iter, pos;
    double *pGT_B, *pD, *pInitProj, *pProj;

    // create pointers to distance and GT_B matrices
    // Setting matrices with values of zero
    d = (double *) malloc((rows * rows) * sizeof(double));
    GT_B = (double *) malloc((rows * rows) * sizeof(double));

    // Matrices of distances
    Tools::EucDist(pData, rows, cols, d);

    // Iterational process
    for (iter = 1; iter <= iterations; iter++) {
        // Creating Guttman tranformation matrix GT_B
        pos = -rows;
        for (i = 0; i < rows; i++) {
            theSum = 0.0;
            pos = pos + rows;
            pGT_B = GT_B + pos;

            // Elements before diagonal
            for (j = 0; j < i; j++, pGT_B++) {
                theSum -= *(pGT_B);
            }
            pGT_B++;

            pD = d + pos + i + 1;

            // Elements after diagonal
            for (j = i + 1; j < rows; j++, pGT_B++, pD++) {
                dist = 0.0;
                for (t = 0; t < projDim; t++) {
                    // Best performance versus pointes or dist = Y * Y
                    Y = *(initProj + projDim * i + t) - *(initProj + projDim * j + t);
                    dist += Y * Y;
                }
                dist = sqrt(dist);

                if (dist > 0.000001) {
                    *(pGT_B) = -*(pD) / dist;
                } else {
                    *(pGT_B) = 0.0;
                }

                *(GT_B + rows * j + i) = *(pGT_B);
                theSum -= *(pGT_B);
            }
            *(GT_B + pos+ i) = theSum;
        }

        // Calculating new projection of points
        // proj = 1 / n * GT_B * proj';
        pos = -rows;
        pProj = proj;
        for (i = 0; i < rows; i++){
            pos = pos + rows;
            // projection coords
            for (t = 0; t < projDim; t++, pProj++) {
                theSum = 0.0;
                pInitProj = initProj + t;
                for (j = 0; (j < rows); j++) {
                    theSum += *(GT_B + pos + j) * *(pInitProj);
                    pInitProj += projDim;
                }
                *(pProj) = theSum / rows;
            }
        }

        // Set proj equal to initProj
        for (i = 0; i < rows * projDim; i++) {
            // projection coords
            *(initProj + i) = *(proj + i);
        }
    }
    free(d);
    free(GT_B);
}

/*
 * 3 times slower than regular Smacof but can achieve similar error (but it can be worse also)
 *
 * May it can be updated to calculate distance matrix after iteration, but points projections after changing in each iteration
 * 
 */
void Mds::SmacofZ(double *pData, int rows, int cols, int iterations, int projDim, double *initProj, double *proj)
{
    double *d, *GT_B, theSum, dist, Y;
    int i, j, t, iter, pos, *index, row;
    double *pGT_B, *pGT_BJ, *pD, *pInitProj, *pProj;
    /* create pointers to distance and GT_B matrices */
    /* Setting matrices with values of zero */
    d = (double *) calloc((rows * rows), sizeof(double));
    GT_B = (double *) calloc((rows * rows), sizeof(double));
    index = (int *) calloc(rows, sizeof(int));

    /* Matrices of distances */
    Tools::EucDist(pData, rows, cols, d);

    /* Creating Guttman tranformation matrix GT_B*/
    pos = -rows;
    for (i = 0; i < rows; i++) {
        theSum = 0.0;
        pos = pos + rows;
        pGT_B = GT_B + pos;

        /* Elements before diagonal */
        for (j = 0; j < i; j++, pGT_B++) {
            theSum -= *(pGT_B);
        }
        pGT_B++;

        pD = d + pos + i + 1;
        /* Elements after diagonal */
        for (j = i + 1; j < rows; j++, pGT_B++, pD++) {
            dist = 0.0;
            for (t = 0; t < projDim; t++) {
                // Best performance versus pointes or dist = Y * Y
                Y = *(initProj + projDim * i + t) - *(initProj + projDim * j + t);
                dist += Y * Y;
            }
            dist = sqrt(dist);

            if (dist > 0.000001) {
                *(pGT_B) = -*(pD) / dist;
            } else {
                *(pGT_B) = 0.0;
            }
            *(GT_B + rows * j + i) = *(pGT_B);
            theSum -= *(pGT_B);
        }
        *(GT_B + pos+ i) = theSum;
    }

    /* Set proj equal to initProj */
    for (i = 0; i < rows * projDim; i++) {
        /* projection coords */
        *(proj + i) = *(initProj + i);
    }

    /* Iterational process */
    for (iter = 1; iter <= iterations; iter++) {
    	/* Calculating new projection of points */
	// proj = 1 / n * GT_B * proj';
        Tools::Shuffle(index, rows, this->GetSeed());

        for (row = 0; row < rows; row++) {
            i = *(index + row);
            pos = i * rows;
            pProj = proj + i * projDim;
            /* projection coords */
            for (t = 0; t < projDim; t++, pProj++) {
                theSum = 0.0;
                pInitProj = proj + t;
                for (j = 0; (j < rows); j++) {
                    theSum += *(GT_B + pos + j) * *(pInitProj);
                    pInitProj += projDim;
		}
                *(pProj) = theSum / rows;
            }

            /* Recalculation of GT_B */
            pGT_B = GT_B + pos;
            pD = d + pos;
            pGT_BJ = GT_B;
            theSum = 0.0;
            for (j = 0; (j < rows); j++, pD++, pGT_B++) {
                dist = 0.0;
                for (t = 0; t < projDim; t++) {
                    // Best performance versus pointers or dist = Y * Y
                    Y = *(proj + projDim * i + t) - *(proj + projDim * j + t);
                    dist += Y * Y;
                }
                dist = sqrt(dist);

                *(pGT_BJ) += *(pGT_B); //
                if (dist > 0.000001) {
                    //*(pGT_B) = -*(pD) / dist;
                    Y = -*(pD) / dist;
                    *(pGT_B) = Y;
                } else {
                    *(pGT_B) = 0.0;
                    Y = 0.0;
                }

                theSum -= Y;
                *(pGT_BJ) -= Y;
                *(pGT_BJ - j + i) = Y;
                pGT_BJ += rows + 1;
            }
            *(GT_B + pos+ i) = theSum;
	}
    }
    free(d);
    free(GT_B);
    free(index);
}
