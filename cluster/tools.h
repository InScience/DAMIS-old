#ifndef TOOLS_H
#define TOOLS_H

#include <fstream> //ifstream
#include <math.h> // sqrt
#include <time.h>
#include "libs/alglib/ap.h"
#include "libs/alglib/correlation.h"

#ifndef DBL_MAX
#define DBL_MAX 1E+300
#endif

class Tools
{
    public:
        static int EucDist(double* pData, int rows, int cols, double* pDataDist);
        static int EucDist(double* pData, int rows, int cols, ap::real_1d_array& pDataDist);
        static int EucDist(double* pData, int rows, int cols, int K, double* pDataDist);
        static void EucDist2Matrix(double* x, double* y, int mx, int my, int n, double* d);
        static void RelativeMdsError(double* d, double* dn, double* proj, double* x, int oldRows, int newRows, int cols, double* fun, double* gr);
        static double MdsError(double* pD, double* pDProj, int iRows);
        static double MinimalWiringError2(double* dataArray, int rows, int cols, double* planeArray, int planeDim, double* wiringArray);
        static double MinimalWiringError(double* dataArray, int rows, int cols, double* planeArray, int planeDim);
        static double BadClassifiedVectors(double* dataArray, int rows, int cols, double* planeArray, int planeDim);
        static double SpearmanCoefficient(double* dataArray, int rows, int cols, double* planeArray, int planeDim);
        static double ApSpearmanCoefficient(double* dataArray, int rows, int cols, double* planeArray, int planeDim);
        static void BubbleSort(double* dataArray, int elements);
        static void BubbleSortDesc(double* dataArray, int elements, int* index);
        static void QuickSortSwap(double *x, double *y);
        static int QuickSortChoosePivot(int i, int j);
        static void QuickSort(double* list, int m, int n);

        static void TimerStart(time_t* tv, int type = 0);
        static double TimerStop(time_t tv, int type = 0);
        static void CalcColsVariance(double* pData, int rows, int cols, double* pVariance, int type = 1);
        static int Shuffle(int* index, int rows = -1, int sead = -1);
    private:

};

#endif
