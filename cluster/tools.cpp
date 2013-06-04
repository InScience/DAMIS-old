#include "tools.h"
//#include <iostream> // laikinai

// Normuota paklaida
double Tools::MdsError(double *pD, double *pDProj, int iRows)
{
    int i, j;
    double theSum, Y, c;
    double *XI, *XJ;

    XI = pD;
    XJ = pDProj;
    c = 0.0;
    theSum = 0.0;
    for (i=0; i<iRows; i++) {
        XI = XI + i + 1;
        XJ = XJ + i + 1;
        for (j=i+1; j<iRows; j++, XI++,XJ++) {
            /* XI = x + i*n; XJ = x + j*n; */
            Y = (*XI)-(*XJ);
            theSum += Y * Y;
            c += (*XI) * (*XI);
        }
    }
    return theSum / c;
}

// Apskaiciuoja euklidinius atstumus tarp vektoriu
// pData - duomenys (vektoriai)
// rows - eilutes
// cols - stulpeliai
// pDist - atstumai
// ----------------
// grazina 0 - nepavyko, 1 - pavyko
int Tools::EucDist(double *pData, int rows, int cols, double *pDataDist)
{
    /*
     d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
    */
    if (rows > 0 && cols > 0) {
        int i, j, k, pos;
        double theSum, Y;
        double *XI, *XJ, *XI0;

        XI = pData;
        pos = 0;
        for (i=0; i<rows; i++,pos=pos+rows) {
            XI0 = XI;
            XJ = XI + cols;
            //*(pDataDist + rows * i + i) = 0;
            *(pDataDist + pos + i) = 0;
            for (j=i+1; j<rows; j++) {
                /* XI = x + i*cols; XJ = x + j*cols; */
                XI = XI0;
                theSum = 0;
                for (k=0; k<cols; k++,XI++,XJ++) {
                    Y = (*XI)-(*XJ);
                    theSum += Y*Y;
                }
                /* *(d++) = sqrt(theSum); */
                //*(pDataDist + rows * i + j) = sqrt(theSum);
                //*(pDataDist + rows * j + i) = *(pDataDist + rows * i + j);
                *(pDataDist + pos + j) = sqrt(theSum);
                *(pDataDist + rows * j + i) = *(pDataDist + pos + j);
            }
        }
        return 1;
    }
    return 0;
}

/**
 * Calculate distacne array
 *
 * @param pData
 * @param rows pData rows
 * @param cols pData cols
 * @param pDataDist size of rows * (rows - 1) / 2 is result
 *
 * @return 1 if success and 0 if not
 */
int Tools::EucDist(double *pData, int rows, int cols, ap::real_1d_array& pDataDist)
{
    /*
     d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
    */
    if (rows > 0 && cols > 0) {
        int i, j, k, pos;
        double theSum, Y;
        double *XI, *XJ, *XI0;

        XI = pData;
        pos = 0;
        for (i=0; i<rows; i++) {
            XI0 = XI;
            XJ = XI + cols;
            //*(pDataDist + rows * i + i) = 0;
            //pDataDist(pos) = 0;
            //pos++;
            for (j=i+1; j<rows; j++, pos++) {
                /* XI = x + i*cols; XJ = x + j*cols; */
                XI = XI0;
                theSum = 0;
                for (k=0; k<cols; k++,XI++,XJ++) {
                    Y = (*XI)-(*XJ);
                    theSum += Y*Y;
                }
                pDataDist(pos) = sqrt(theSum);
            }
        }
        return 1;
    }
    return 0;
}


/* Euclidean distance between two matrices for DMA  */
int Tools::EucDist(double *pData, int rows, int cols, int K, double *pDataDist)
{
    /*
     d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
     */
    int i, j, k, pos;
    double theSum, Y;
    double *XI, *XJ, *XI0;

    XI = pData;
    pos = 0;
    for (i=0; i<rows; i++,pos=pos+rows) {
        XI0 = XI;
        XJ = XI+cols;
        *(pDataDist + pos+ i) = 0.0;
        for (j=i+1; (j<rows && j <= i + K); j++) {
        /* XI = x + i*cols; XJ = x + j*cols; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<cols;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y;
            }
            *(pDataDist + pos + j) = sqrt(theSum);
            *(pDataDist + rows * j + i) = *(pDataDist + pos + j);
        }
    }
}

/** Uzfiksuoja laika
 *
 * clock yra long int tipo ir todėl esant dideliems CLOCKS_PER_SEC skaitliukas
 * persisuka ir skaičiuoja per naują
 * 
 * @param tv - laikass
*/
void Tools::TimerStart(time_t *tv, int type)
{
    if (type == 0) {
        (*tv) = clock();
    } else {
        (*tv) = time(NULL);
    }
}

// Suskaiciuoja laiko pokyti
// tv - laikas
// --------------------------
// Grazina laiko skirtuma tarp tv ir dabar
double Tools::TimerStop(time_t tv, int type)
{
    if (type == 0) {
        // Po 35min laikas tampa neigimas
        // todėl prie gauto laiko reikia pridėti 2147,483647 s
        return (double)((0.0 + clock() - tv) / CLOCKS_PER_SEC);
    } else {
        // laikas sekundžių tikslumu
        time_t end = time(NULL);
        return difftime(end, tv);
    }
}

/**
 * Apskaiciuoja eiluciu dispersijas
 * Naive algorithm
 *
 * @param pData - data matrix
 * @param rows  - a number of rows
 * @param cols  - a number of cols
 * @param pVar  - array containing cols variance (/n)
 * @param type  - 1 naive algorithm, 2 - two path algorithm, 3 - compensated algorithm
 */
void Tools::CalcColsVariance(double *pData, int rows, int cols, double* pVar, int type)
{
    int i, j;
    double *X, *Y, *ave, tmp;

    if (rows < 1 ||cols < 1) {
        return;
    }

    //ave = (double *) malloc (cols * sizeof(double));
    ave = (double *) calloc (cols, sizeof(double)); // Array inital values is 0
    for (i = 0; i < cols; i++) {
        //*(ave + i) = 0.0;
        *(pVar + i) = 0.0;
    }

    switch (type) {
        default:
        // Naive algorithm
        case 1: // Fastest
            // Mean and X^2
            X = pData;
            j = 0;
            for (i = 0; i < cols * rows; i++, X++, j++) {
                if (j >= cols) {
                    j = 0;
                }
                *(ave + j)  += *(X);
                *(pVar + j) += *(X) * *(X);
            }

            // Variance
            X = ave;
            Y = pVar;
            for (i = 0; i < cols; i++, X++, Y++) {
                *(X) = *(X) / rows;
                *(Y) = *(Y) / rows - *(X) * *(X);
            }
        break;
        // Two path algorithm
        case 2: // 1,42 times slower than 1
            // Mean
            X = pData;
            j = 0;
            for (i = 0; i < cols * rows; i++, X++, j++) {
                if (j >= cols) {
                    j = 0;
                }
                *(ave + j)  += *(X);
            }

            X = ave;
            for (i = 0; i < cols; i++, X++) {
                *(X) = *(X) / rows;
            }

            // Variance
            X = pData;
            j = 0;
            for (i = 0; i < cols * rows; i++, X++, j++) {
                if (j >= cols) {
                    j = 0;
                }
                tmp = (*(X) - *(ave + j));
                *(pVar + j) += tmp * tmp;
            }

            X = pVar;
            for (i = 0; i < cols; i++, X++) {
                *(X) = *(X) / rows;
            }
        break;
        // Compensated
        case 3: // 1,78 times slower than 1
            double *sumc;
            sumc = (double *) malloc (cols * sizeof(double));
            // Mean
            X = pData;
            j = 0;
            for (i = 0; i < cols * rows; i++, X++, j++) {
                if (j >= cols) {
                    j = 0;
                }
                *(ave + j) += *(X);
            }

            X = ave;
            for (i = 0; i < cols; i++, X++) {
                *(X) = *(X) / rows;
                *(sumc + i) = 0.0;
            }

            // Variance
            X = pData;
            j = 0;
            for (i = 0; i < cols * rows; i++, X++, j++) {
                if (j >= cols) {
                    j = 0;
                }
                tmp = (*(X) - *(ave + j));
                *(pVar + j) += tmp * tmp;
                *(sumc + j) += tmp;
            }

            X = pVar;
            Y = sumc;
            for (i = 0; i < cols; i++, X++, Y++) {
                *(X) = (*(X) - (*(Y) * *(Y)) / rows) / rows;
            }
            free(sumc);
        break;
    }
    free(ave);
}

/**
 *  sumaisome index masyvą
 */
int Tools::Shuffle(int *index, int rows, int sead)
{
    int i, tmp, index1, index2;

    for (i = 0; i<rows; i++) {
        *(index + i) = i;
    }

    if (sead < 1) {
        srand((unsigned)time( NULL ));
    } else {
        srand(sead);
    }

    for(i = 0; i <= 5 * rows; i++)
    {
        index1 = (int) fmod(rand(), rows + 0.0);
        index2 = (int) fmod(rand(), rows + 0.0);

        tmp =  *(index + index2);
        *(index + index2) = *(index + index1);
        *(index + index1) = tmp;
    }

    return 1;
}

/* Euclidean distance between two matrices */
void Tools::EucDist2Matrix(double *x, double *y, int mx, int my, int n, double *d)
{
    /*
     d = sqrt(sum((XI-XJ).^2, 2));            % Euclidean
     */
    int i, j, k;
    double theSum, Y;
    double *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<mx; i++) {
        XI0 = XI;  // Dideja vis per viena vektoriu
        XJ = y;    // Grazinamas prie pirmo Y matricos vektoriaus
        for (j=0; j<my; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y;
            }
            *(d + mx * j + i) = sqrt(theSum);  // Grazina atsumu matrica [mx x my]
        }
    }
}

void Tools::RelativeMdsError(double *d, double *dn, double *proj, double *x, int oldRows, int newRows, int cols, double *fun, double *gr)
{
    int i, j, k, t;
    double *projD, *xprojD, *inxI, *inxJ, Y;
    double sumNewNew, sumNewOld;

    projD = (double *) malloc(oldRows * newRows * sizeof(double));

    /* Calculation of sumNewNew.It's error between new point projection distances and in original data */
    sumNewNew = 0.0; // Pirmas elementas 0
    if (newRows > 1) {
        xprojD = (double *) malloc(newRows * newRows *  sizeof(double));
        Tools::EucDist(x, newRows, cols, xprojD);
        inxI = dn;  // To speed up
        inxJ = xprojD;
        for (i=0; i<newRows; i++) {
            for (j=0; j<newRows; j++, inxI++, inxJ++) {
                //sumNewNew += pow((*(dn + i * newRows + j) - *(xprojD + i * newRows + j)), 2) ;
                sumNewNew += pow((*(inxI) - *(inxJ)), 2);
                t = i * newRows + j;
                for (k=0; k<cols; k++) {
                     if (*(xprojD + t) > 0.000001)
                        *(gr + k * newRows + j) += 2 * (*(dn + t) - *(xprojD + t)) / *(xprojD + t) * (*(x + i * cols + k) - *(x + j * cols + k));
                }
            }
        }
        sumNewNew /= 2;
    }
    
    /* Calculation of sumNewOld. It's error between new and old points projections distances and in original data */
    sumNewOld = 0.0;
    EucDist2Matrix(proj, x, oldRows, newRows, cols, projD);
    
    inxI = d;  // To speed up
    inxJ = projD;
    for (i=0; i<oldRows; i++) {
        for (j=0; j<newRows; j++, inxI++, inxJ++) {
            //sumNewOld += pow((*(d + i * newRows + j) - *(projD + i * newRows + j)), 2) ;
            //sumNewOld += pow((*(inxI) - *(inxJ)), 2) ;
            Y = *(inxI) - *(inxJ);
            sumNewOld += Y * Y;
            
            /* Gradient */
            t = j * oldRows + i;
            if (*(projD + t) > 0.000001) {
                for (k=0; k<cols; k++) {
                    *(gr + k * newRows + j) += 2 * (*(d + t) - *(projD + t)) / *(projD + t) * (*(proj + i * cols + k) - *(x + j * cols + k));
                }
            }
        }
    }

    *(fun) = sumNewNew + sumNewOld ;
}

/**
 *
 * Minimalaus jungimo paklaida
 *
 * Grazina klaida ir sujungimu matrica
 */
double Tools::MinimalWiringError2(double *dataArray, int rows, int cols, double *planeArray, int planeDim, double *wiringArray)
{
    double *d, *d2;
    int i, j;

    // Iskiriame masyvams vieta
    // Daugiamaciams masyvams
    d = (double *) calloc (rows * rows, sizeof(double));

    d2 = (double *) calloc (rows * rows, sizeof(double));

    Tools::EucDist(dataArray, rows, cols, d);

    Tools::EucDist(planeArray, rows, planeDim, d2);

    // Prim'o algoritmas minimalaus jungimo medziui rasti
    double *key;
    double error;
    int *p;
    int u, v, pos;

    // Iskiriame masyvams vieta
    key = (double *) calloc (rows, sizeof(double));
    p = (int *) calloc (rows, sizeof(int));

    // Priskiriame pradines reiksmes
    for (i = 0; i < rows; i++) {
        key[i] = DBL_MAX;
        p[i] = -1;
    }

    // Pasirenkame elementa no kurio pradesime
    u = 0;
    key[u] = 0.0;
    p[u] = -1;

    pos = 0;
    // Reikia kad paimsime visas krastines
    for (u = 0, pos = 0; u < rows - 1; u++, pos = pos + rows) {
        //u = i;
        // Pradedame ieskoti atstumo ir kaimynu
        for (v = u + 1; v < rows; v++) {
            if (*(d + pos + v) < *(key + v)) {
                *(p + v) = u;
                *(key + v) = *(d + pos + v);
            }
        }
    }

    error = 0.0;
    // Konstruojame rezultata
    for (i = 0, pos = 0; i < rows - 1; i++, pos = pos + rows) {
        for (j = i + 1; j < rows; j++) {
            if (*(p + j) == i) {
                *(wiringArray + pos + j) = 1.0;
                *(wiringArray + j * rows + i) = 1.0;
                error += *(d2 + pos + j);
            }
        }
    }

    free(d);
    free(d2);
    free(key);
    free(p);

    return error;
}

/**
 *
 * Minimalaus jungimo paklaida
 *
 * Grazina klaida
 */
double Tools::MinimalWiringError(double *dataArray, int rows, int cols, double *planeArray, int planeDim = 2)
{
    double *d, *d2;
    int i, j;

    // Iskiriame masyvams vieta
    // Daugiamaciams masyvams
    d = (double *) calloc (rows * rows, sizeof(double));

    d2 = (double *) calloc (rows * rows, sizeof(double));

    Tools::EucDist(dataArray, rows, cols, d);

    Tools::EucDist(planeArray, rows, planeDim, d2);

    // Prim'o algoritmas minimalaus jungimo medziui rasti
    double *key;
    double error;
    int *p;
    int u, v, pos;

    // Iskiriame masyvams vieta
    key = (double *) calloc (rows, sizeof(double));
    p = (int *) calloc (rows, sizeof(int));

    // Priskiriame pradines reiksmes
    for (i = 0; i < rows; i++) {
        key[i] = DBL_MAX;
        p[i] = -1;
    }

    // Pasirenkame elementa no kurio pradesime
    u = 0;
    *(key + u) = 0.0;
    *(p + u) = -1;

    // Reikia kad paimsime visas krastines
    for (u = 0, pos = 0; u < rows - 1; u++, pos = pos + rows) {
        // Pradedame ieskoti atstumo ir kaimynu
        for (v = u + 1; v < rows; v++) {
            if (*(d + pos + v) < *(key + v)) {
                *(p + v) = u;
                *(key + v) = *(d + pos + v);
            }
        }
    }

    error = 0.0;
    // Konstruojame rezultata
    for (i = 0, pos = 0; i < rows - 1; i++, pos = pos + rows) {
        for (j = i + 1; j < rows; j++) {
            if (*(p + j) == i) {
                error += *(d2 + pos + j);
            }
        }
    }
    free(d);
    free(d2);
    free(key);
    free(p);

    return error;
}

// Neislaikyti atstumai
double Tools::BadClassifiedVectors(double *dataArray, int rows, int cols, double *planeArray, int planeDim = 2)
{
    double *d, *d2;
    int i, i2, j, j2, k, k2;
    //double num, total;

    // Iskiriame masyvams vieta
    // Daugiamaciams masyvams
    d = (double *) calloc (rows * rows, sizeof(double));

    d2 = (double *) calloc (rows * rows, sizeof(double));

    Tools::EucDist(dataArray, rows, cols, d);

    Tools::EucDist(planeArray, rows, planeDim, d2);


    // Surasome elementus i viena masyva
    double total, num, t;

    t = 0.5 * rows * (rows - 1);
    total = ((t * (t - 1) / 2) - 1) * (t * (t - 1) / 2) / 2;

    // Nustatome atstumus tarp vektoriu - n maciu ir dvimaciu
    // n^4 sudetingumas skaiciuoja labai letai
    for (i = 0, k = 0; i < rows - 1; i++, k = k + rows) {
        for (j = i + 1; j < rows; j++) {
            for (i2 = 0, k2 = 0; i2 < rows - 1; i2++, k2 = k2 + rows) {
                for (j2 = i2 + 1; j2 < rows; j2++) {
                    if (*(d + k + j) < *(d + k2 + j2)) {
                        if (*(d2 + k + j) <= *(d2 + k2 + j2)) {
                            num += 1.0;
                        }
                    }
                }
            }
        }
    }

    // Isvalome atminti nuo dinaminiu kintmųjų
    free(d);
    free(d2);
    
    return (num / total);
}

/**
 *
 * Speramano koeficientas
 *
 * Dubliuoja ApSpearmanCoefficient, nes senesnis algoritmas veike ilgai ir nekorektiskai
 * 
 */
double Tools::SpearmanCoefficient(double *dataArray, int rows, int cols, double *planeArray, int planeDim = 2)
{ 
    return Tools::ApSpearmanCoefficient(dataArray, rows, cols, planeArray, planeDim);
}

/**
 * procedure bubbleSort( A : list of sortable items ) defined as:
  n := length( A )
  do
    swapped := false
    for each i in 0 to n - 1  inclusive do:
      if A[ i ] > A[ i + 1 ] then
        swap( A[ i ], A[ i + 1 ] )
        swapped := true
      end if
    end for
    n := n - 1
  while swapped
end procedure
*/
void Tools::BubbleSort(double* dataArray, int elements)
{
    bool swapped;
    double tmp;

    do {
        swapped = false;
        for (int i = 1; i < elements; i++) {
            if (dataArray[i] > dataArray[i-1]) {
                tmp = dataArray[i-1];
                dataArray[i-1] = dataArray[i];
                dataArray[i] = tmp;
                swapped = true;
            }
        }
        elements--;
    } while (swapped);
}

/**
 * Make bubble sorting of dataArray and return index matrix which contain
 * information about index changes. Data is in descending order
 *
 * @param double 1 dimensional real numbers array dataArray
 * @param int elements how many
 * @param integers one dimensional array it shut be inicialized before
 *
 * @return void
 */
void Tools::BubbleSortDesc(double* dataArray, int elements, int* index)
{
    bool swapped;
    double tmp;
    int t;

    // Set index array values
    for (int i = 0; i<elements; i++) {
        *(index + i) = i;
    }

    do {
        swapped = false;
        for (int i = 1; i < elements; i++) {
            if (dataArray[i] > dataArray[i-1]) {
                tmp = dataArray[i-1];
                dataArray[i-1] = dataArray[i];
                dataArray[i] = tmp;
                // Change indexes
                t = index[i-1];
                index[i-1] = index[i];
                index[i] = t;
                swapped = true;
            }
        }
        elements--;
    } while (swapped);
}

/**
 *
 * Speramano koeficientas
 *
 */
double Tools::ApSpearmanCoefficient(double *dataArray, int rows, int cols, double *planeArray, int planeDim = 2)
{
    double coef = 0.0;
    int size;
    ap::real_1d_array d_;
    ap::real_1d_array d2_;

    size = (int ) rows * (rows - 1) / 2;

    // Surasome elementus i viena masyva
    d_.setlength(size);
    d2_.setlength(size);

    // Apskaičiuojame Euklido
    Tools::EucDist(dataArray, rows, cols, d_);
    Tools::EucDist(planeArray, rows, planeDim, d2_);

    // apskaiciuojame korelacija naudodamiesi ALGLIB
    //if (rows < 4000) {
        coef = spearmanrankcorrelation(d_, d2_, size);
   // } else {
   //     coef = -1;
   //     coef = pearsoncorrelation(d_, d2_, size);
   // }

    return coef;
}


void Tools::QuickSortSwap(double *x, double *y)
{
   double temp;
   temp = *x;
   *x = *y;
   *y = temp;
}

int Tools::QuickSortChoosePivot(int i, int j)
{
   return int((i+j) / 2);
}

void Tools::QuickSort(double* list, int m, int n)
{
    double key;
    int i, j, k;
    if (m < n) {
        k = Tools::QuickSortChoosePivot(m, n);
        Tools::QuickSortSwap(&list[m], &list[k]);
        key = list[m];
        i = m + 1;
        j = n;
        while (i <= j) {
            while ((i <= n) && (list[i] <= key))
                i++;
            while ((j >= m) && (list[j] > key))
                j--;
            if (i < j)
                Tools::QuickSortSwap(&list[i], &list[j]);
        }
      // swap two elements
      Tools::QuickSortSwap(&list[m], &list[j]);
      // recursively sort the lesser list
      Tools::QuickSort(list, m, j-1);
      Tools::QuickSort(list, j+1, n);
   }
}
