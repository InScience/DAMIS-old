///////////////////////////////////////////////////////////
//  PCA.cpp
//  Implementation of the Class PCA
//  Created on:      07-Lie-2013 20:07:31
//  Original author: Povilas
///////////////////////////////////////////////////////////

/*! \class PCA
    \brief A class of methods and attributes for PCA algorithm.
 */

#include <math.h>
#include "PCA.h"
#include "Statistics.h"
#include "DistanceMetrics.h"
#include "alglib/dataanalysis.h"

PCA_::PCA_(){

}

PCA_::~PCA_(){

}

PCA_::PCA_(double d, bool disp){
    if (disp)
    {
        setProjectionDimension((int)d);
        initializeProjectionMatrix();
        ProjectXMatrix();
    }
    else
    {
        int n = X.getObjectCount();
        int m = X.getObjectAt(0).getFeatureCount();
        int dd = 0;
        PCA_ pca(X, m);
        ObjectMatrix Y_visi = pca.getProjection();
        this->eigValues = pca.getEigenValues();
        double wholeSum = 0.0, tempSum = 0.0;

        for (int i = 0; i < m; i++)
            wholeSum += eigValues(i);

        for (int i = 0; i < m; i++)
        {
                tempSum += eigValues(dd);
                dd++;
            if (double(tempSum / wholeSum)  > (double) d / 100.0)
                break;
        }

        setProjectionDimension(dd);
        initializeProjectionMatrix();

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < dd; j++)
                Y.updateDataObject(i, j, Y_visi.getObjectAt(i).getFeatureAt(j));
        }
    }
}

PCA_::PCA_(ObjectMatrix Xmatrix, double disp){
    X = Xmatrix;
    dispPart = disp;
    int n = X.getObjectCount();
    int m = X.getObjectAt(0).getFeatureCount();
    int dd = 0;
    PCA_ pca(X, m);
    ObjectMatrix Y_visi = pca.getProjection();
    this->eigValues = pca.getEigenValues();
    double wholeSum = 0.0, tempSum = 0.0;

    for (int i = 0; i < m; i++)
        wholeSum += eigValues(i);

    for (int i = 0; i < m; i++)
    {
        tempSum += eigValues(i);
        dd++;
        if ((tempSum / wholeSum) > disp)
            break;
    }
    d = dd;
    setProjectionDimension(d);
    initializeProjectionMatrix();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < d; j++)
            Y.updateDataObject(i, j, Y_visi.getObjectAt(i).getFeatureAt(j));
    }
}

PCA_::PCA_(ObjectMatrix objMatrix, int dim){
    setProjectionDimension(dim);
    X = objMatrix;
    initializeProjectionMatrix();
    ProjectXMatrix();
}

double PCA_::getStress()
{
    //double stress = 0.0;
    double stress = 0.0;
    int m = X.getObjectCount();
    double distX = 0.0;
    double distY = 0.0;

    for (int i = 0; i < m - 1; i++)
    {
        for (int j = i + 1; j < m; j++)
        {
            distX = DistanceMetrics::getDistance(X.getObjectAt(i), X.getObjectAt(j), EUCLIDEAN);
            distY = DistanceMetrics::getDistance(Y.getObjectAt(i), Y.getObjectAt(j), EUCLIDEAN);
            stress += pow(distX - distY, 2);
        }
    }

    return stress * DimReductionMethod::getStressWeight();

    //return stress;
}

ObjectMatrix PCA_::getProjection(){

    return Y;
}

void PCA_::toDataType(){
    int m = X.getObjectCount();
    int n = X.getObjectAt(0).getFeatureCount();
    alglibX.setlength(m, n);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            alglibX(i,j) = X.getObjectAt(i).getFeatureAt(j);
}

void PCA_::ProjectXMatrix()
{
    PCA_::toDataType();
    int m = X.getObjectCount();
    int n = X.getObjectAt(0).getFeatureCount();
    std::vector<double> X_vid;
    X_vid.reserve(0);
    //double X_vid[n],
    double tmp = 0.0;
    double wholeDisp = 0.0, tarpDisp = 0.0;

    for (int i = 0; i < n; i++)
        X_vid.push_back(Statistics::getAverage(X, i)) ;

    alglib::ae_int_t info;
    alglib::real_2d_array eigVectors;
 //   alglibX.setcontent(m,n,)
    pcabuildbasis(alglibX, m, n, info, eigValues, eigVectors);

    if (info == 1)
    {
        for (int i = 0; i < m; i++) // objektu skaicius
            for (int j = 0; j < d; j++) //i kuria erdve projekuojam
            {
                tmp = 0.0;
                for (int k = 0; k < n; k++) // objekto atriutu skaicius
                    tmp += (alglibX(i,k) - X_vid.at(k)) * eigVectors[k][j];
                Y.updateDataObject(i, j, tmp);
            }

        for (int i = 0; i < d; i++)
            tarpDisp += eigValues[i];

        for (int i = 0; i < n; i++)
            wholeDisp += eigValues[i];

        dispPart = tarpDisp / wholeDisp;
    }
}

void PCA_::fromDataType(){
    int m = X.getObjectCount();
    int n = X.getObjectAt(0).getFeatureCount();

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            X.updateDataObject(i, j, alglibX(i,j));

}

int PCA_::getDimension()
{
    return d;
}

double PCA_::getDispersionPart()
{
    return dispPart;
}

alglib::real_1d_array PCA_::getEigenValues()
{
    return eigValues;
}
