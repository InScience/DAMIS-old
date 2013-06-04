/*************************************************************************
Copyright (c) 2008, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#include "alglib_pca.h"

/*************************************************************************
Principal components analysis

Subroutine  builds  orthogonal  basis  where  first  axis  corresponds  to
direction with maximum variance, second axis maximizes variance in subspace
orthogonal to first axis and so on.

It should be noted that, unlike LDA, PCA does not use class labels.

INPUT PARAMETERS:
    X           -   dataset, array[0..NPoints-1,0..NVars-1].
                    matrix contains ONLY INDEPENDENT VARIABLES.
    NPoints     -   dataset size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1

�������� ���������:
    Info        -   return code:
                    * -4, if SVD subroutine haven't converged
                    * -1, if wrong parameters has been passed (NPoints<0,
                          NVars<1)
                    *  1, if task is solved
    S2          -   array[0..NVars-1]. variance values corresponding
                    to basis vectors.
    V           -   array[0..NPoints-1,0..NVars-1]
                    matrix, whose columns store basis vectors.

*************************************************************************/
void calculatePcaProjection(const ap::real_2d_array& x,
    int npoints,
    int nvars,
    int nprojdim,
    double *proj)
{
    ap::real_2d_array apProj;
    int info;
    ap::real_1d_array s2, work; // it's unclear why where "work" is used for
    ap::real_2d_array v;

    ap::ap_error::make_assertion(npoints>=0, "npoints shut be more than 0!");
    ap::ap_error::make_assertion(nvars>=0, "nvars shut be more than 0!");
    ap::ap_error::make_assertion(nprojdim>=0&&nprojdim<=nvars, "dim shut be more than 0 and less than number of vnars!");

    // Calculate svd decomposition of dataset and get V^T = v
    pcabuildbasis(x, npoints, nvars, info, s2, v);

    // multiplly X * v
    int i = ap::maxint(npoints, nvars);
    work.setbounds(1, i);
    // Kai skaiciuojame visa matrica
    //apProj.setbounds(0, npoints-1, 0, nvars-1);
    //matrixmatrixmultiply(x, 0, npoints-1, 0, nvars-1, false, v, 0, nvars-1, 0, nvars-1, false, 1, apProj, 0, npoints-1, 0, nvars-1, 0, work);
    // Skaiciauojame tik tiek elementu kiek reikia
    apProj.setbounds(0, npoints-1, 0, nprojdim-1);
    matrixmatrixmultiply(x, 0, npoints-1, 0, nvars-1, false, v, 0, nvars-1, 0, nprojdim-1, false, 1, apProj, 0, npoints-1, 0, nprojdim-1, 0, work);

    // Perverciam i musu formata
    // proj turi buti inicializuotas pradzioje pries paduodant i funkcija
    for (int i=0; i<npoints; i++) {
        for (int j=0; j<nprojdim; j++) {
            *(proj + i * nprojdim +j) = apProj(i, j);
        }
    }
}




