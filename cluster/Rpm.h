/* 
 * File:   Rpm.h
 * Author: Virginijus
 *
 * Created on Ketvirtadienis, 2010, Spalio 28, 15.11
 */
#ifndef RPM_H
#define	RPM_H

#include "Projection_Abstract.h"
#include <stddef.h> //NULL
#include <cstdlib> //relloc, malloc

#ifndef INIT_ARRAY_SIZE
#define INIT_ARRAY_SIZE 1000
#endif

#define PI 3.14159265

using namespace std;

class Rpm : public Projection_Abstract {
    //friend int Mds::InitProjectionVectors(double* pData, int method, int rows, int cols, int projDim, double* pPlaneArr);
public:
    Rpm();
    virtual ~Rpm();

    double getError();
    void   train(int);
    void   distances(double*);
    double Hij(int, int, int);
    
    void   movePointInRectangle(double, double);
    double distance(double*, double*);

    static void error();



    
private:
    
    double _error;
    
    double _width;
    double _height;
    
    double _p;

};

#endif	/* RPM_H */

