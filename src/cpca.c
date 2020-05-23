#include <stdlib.h>
#include <stdio.h>
#include "linal.h"
#include "cpca.h"
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
Returns 1 for triangular matrices
*/
char istriangular(matrix* mat)
{
    if(!mat)return 0;
    unsigned long int x,y,h;
    for(y=1;y<mat->h;y++)
    {
        h=mat->w*y;
        for(x=0;x<y;x++)
            if(mat->data[h+x]*mat->data[h+x]>0.0000000000001)return 0;
    }
    return 1;
}

/*
Returns all the eigenpairs for a given matrix <mat>
*/
eigenpairs* qr_eigenpairs_square(matrix* mat)
{
    if(!mat||mat->h!=mat->w)return NULL;
    eigenpairs* ret=malloc(sizeof(eigenpairs));
    ret->eigenvalues=malloc(sizeof(double)*mat->w);
    qr_decomp* QR=qr_decompose(mat);
    matrix* RQ=multiply(QR->R,QR->Q);
    matrix* Q=Matrix(QR->Q->w,QR->Q->h,QR->Q->data),*tmp;
    unsigned long int i;
    while(!istriangular(RQ))
    {
        freeMatrix(&QR->R);
        freeMatrix(&QR->Q);
        free(QR);
        QR=qr_decompose(RQ);
        freeMatrix(&RQ);
        tmp=multiply(Q,QR->Q);
        freeMatrix(&Q);
        Q=tmp;
        RQ=multiply(QR->R,QR->Q);
    }
    ret->eigenvectors=Q;
    for(i=0;i<mat->w;i++)
        ret->eigenvalues[i]=RQ->data[i+i*mat->w];
    freeMatrix(&QR->Q);
    free(QR);
    freeMatrix(&RQ);
    return ret;
}

matrix* square_from_vect(double* vect,unsigned long int len)
{
    if(!vect)return NULL;
    unsigned long int i;
    matrix* ret=Matrix(len,len,NULL);
    for(i=0;i<len;i++)
        ret->data[i+i*len]=vect[i];
    return ret;
}

SVD* qr_svd(matrix* mat)
{
    SVD* ret=NULL;
    if(!mat)return NULL;
    printf("IN SVD\n");
    ret=malloc(sizeof(SVD));
    matrix *trans=transposed(mat),*sq=multiply(trans,mat);
    eigenpairs* eigen=qr_eigenpairs_square(sq);
    ret->U=transposed(eigen->eigenvectors);
    ret->S=square_from_vect(eigen->eigenvalues,sq->w);
    printf("MtM:\n");
    printmat(sq);
    printf("U:\n");
    printmat(ret->U);
    printf("S:\n");
    printmat(ret->S);
    printf("RETURNING FROM SVD\n");
    return ret;
}