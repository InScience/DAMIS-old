/* 
 * File:   Rpm.cpp
 * Author: Virginijus
 * 
 * Created on Ketvirtadienis, 2010, Spalio 28, 15.11
 */
#include "Rpm.h"

Rpm::Rpm() {
    this->_error          = -1.0;
    this->_cols           = 0;
    this->_rows           = 0;
    this->_projDim        = 2;        // fixed
    this->_projInitMethod = 1; // random
    this->_width          = 1.0; // Rectangle width
    this->_height         = 1.0; // Rectangle height
    this->_p              = 0;
            

    this->_pData = (double *) malloc(INIT_ARRAY_SIZE * sizeof(double));
    this->_pProj = (double *) malloc(INIT_ARRAY_SIZE * sizeof(double));
}

Rpm::~Rpm() {
}

void Rpm::train(int iterations) {
    if (this->_pData == NULL && strlen(this->_pFileName) > 0) {
        this->loadData(_pFileName);
    }

    // Start training
    if (this->_rows > 0 && this->_cols > 0) {
        this->_pProj = (double*) realloc(this->_pProj, this->_rows * this->_projDim * sizeof(double));
        // Init projection vectors
        if (this->initProjectionVectors(this->_pData, this->_projInitMethod, this->_rows, this->_cols, this->_projDim, this->_pProj) > 0) {
            double *dataD, *projD;
            // Setting matrice values to zeros
            dataD = (double *) malloc(this->_rows * this->_rows * sizeof(double));
            projD = (double *) malloc(this->_rows * this->_rows * sizeof(double));
            //@TODO speed up
            for (int i=0; i<this->_rows; i++) {
                for (int j=i; j<this->_rows; j++) {
                    *projD(i * this->_rows + j) = this->distance(this->_pProj[i], this->_pProj[j]);
                    *projD(j * this->_rows + i) = *projD(i * this->_rows + j);
                }
            }

            // Matrice of data distances
            Tools::EucDist(this->_pData, this->_rows, this->_cols, dataD);
            
            // Matrice of projetion distances
            this->distances(projD);

            // Energy function first derivate by distance and set to zero
            double *Ep1X = (double *) malloc(this->_rows * sizeof(double));
            double *Ep1Y = (double *) malloc(this->_rows * sizeof(double));
            double *Ep2X = (double *) malloc(this->_rows * sizeof(double));
            double *Ep2Y = (double *) malloc(this->_rows * sizeof(double));

            // Energy function first derivate by distance
            // Fij = -dataD ./ (projD.^(p + 1));
            double *Fij = (double *) malloc(this->_rows * this->_rows * sizeof(double));
            for (int i=0; i<this->_rows; i++) {
                for (int j=i; j<this->_rows; j++) {
                    *Fij(i * this->_rows + j) = -dataD / pow(projD[i * this->_rows + j], this->_p);
                    *Fij(j * this->_rows + i) = *Fij(i * this->_rows + j);
                }
            }

            // Distance first derivate by coordinates
            for (int i=0; i<this->_rows; i++) {
                for (int j=0; j<this->_rows; j++) {
                    //First derivate
                        Ep1X[i] = Ep1X[i] + Hij(i, j, 0) * Fij[i * this->_rows + j];
                        Ep1Y[i] = Ep1Y[i] + Hij(i, j, 1) * Fij[i * this->_rows + j];

                    //Second derivate
                    if (abs(projD(i, j)) > 1e-6) {
                        Ep2X[i] = Ep2X[i] + Fij[i * this->_rows + j] / projD[i * this->_rows + j];
                    }
                }
                Ep2X[i] = -(this->_p + 1) * Ep2X[i];
                Ep2Y[i] = Ep2X[i];
            }
        }
        
        // Training loop
        for (int iter = 0; iter < iterations; iter++) {
            
        }
    } 
}





/**
 * Return Hij value depend on coordinate of projection point
 * 
 * @param int idxI
 * @param int idxJ
 * @param int coordinate
 * @return double
 */
double Rpm::Hij(int idxI, int idxJ, int coordinate) {
// Distance first derivate by coordinates X and Y
//       HijX = sin(2 * pi / w * (xi - xj) ); 
//      HijY = sin(2 * pi / h * (proj_i(2) - proj_j(2)) ); 
    double xi = *(this->_pProj + this->_rows * idxI + coordinate);
    double xj = *(this->_pProj + this->_rows * idxJ + coordinate);
    double w;
    if (coordinate == 0) {
        w = this->_width;
    } else {
        w = this->_height;
    }
    
      
    if ((abs(xi - xj)) < w / 2)  {
        if (xi > xj) {
            return 1 / w;
        } else {
            if (xi < xj) {
                return -1 / w;
            } else {
                return  1 / w;    // if xi = xj
            }
        }
    } else 
        if (abs(xi - xj) > w / 2) {
            if (xi > xj) {
                return -1 / w;
            } else 
                if (xi < xj) {
                    return 1 / w;
                }
        } else {
            return 1 / w;
        }
}


void Rpm::error() {

}



double Rpm::getError() {
    if (this->_error < 0) {
        
    }
    return this->_error;
}

/**
 * This procedure recalculate point coordinates (2dim) to keep him in rectangle [w x h]
 * 
 * @param double x
 * @param double y
 */
void Rpm::movePointInRectangle(double x, double y) {
    double w = this->_width;
    double h = this->_height;
    
    while ((x < 0) || (x >= w)) {
        if (x < 0) {
            x = x + w;
        } else {
            x = x - w;
        }
    }

    while ((y < 0) || (y >= h)) {
        if (y < 0) {
            y = y + h;
        } else {
            y = y - h;
        }
    }
}

/**
 * This function calculate modified distance between voctor1 and vector2.
 * For 2 dimensional vectors only 
 * @param vector_i
 * @param vector_j
 * @return double distance beatween to
 */
double Rpm::distance(double* vector_i, double* vector_j) {
    return this->_width  / (4) * (1 - cos(2 * PI / this->_width  * (vector_i(1) - vector_j(1)))) +
           this->_height / (4) * (1 - cos(2 * PI / this->_height * (vector_i(2) - vector_j(2))));
}

/**
 * Calculate matrix of distances
 * 
 * @param double* projD this matrix is returned
 */
void Rpm::distances(double* projD) {
    for (int i=0; i<this->_rows; i++) {
        for (int j=i; j<this->_rows; j++) {
            *projD(i * this->_rows + j) = this->distance(this->_pProj[i], this->_pProj[j]);
            *projD(j * this->_rows + i) = *projD(i * this->_rows + j);
        }
    }
}