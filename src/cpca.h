#ifndef CPCA_H
#define CPCA_H
#include "linal.h"
/*
A naïve SVD implementation for PCA in BMP images
Copyright (c) 2020 Amélia O. F. da S.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
Eigenpair structure
*/
typedef struct eigenpairs{
    matrix* eigenvectors;
    double* eigenvalues;
}eigenpairs;

/*
Returns 1 for triangular matrices
*/
char istriangular(matrix* mat);

/*
Performs a naïve QR algorithm and returns the eigenpairs for a given square matrix
*/
eigenpairs* qr_eigenpairs_square(matrix* mat);

typedef struct SVD{
    matrix* U;
    matrix* S;
    matrix* VT;
}SVD;

/*
Performs a naïve SVD
*/
SVD* qr_svd(matrix* mat);
#endif