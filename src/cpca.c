#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "linal.h"
#include "cpca.h"

#include "../test/gfx.h"
void showmat(matrix* mat,char d)
{
    if(!mat)return;
    double min=__DBL_MAX__,max=__DBL_MIN__,map,col;
    unsigned int x,y,nx,ny;
    unsigned long int h;
    unsigned char c;
    for(y=0;y<mat->h;y++)
    {
        h=y*mat->w;
        for(x=0;x<mat->w;x++)
        {
            col=mat->data[h+x];
            min=col<min?col:min;
            max=col>max?col:max;
        }
    }
    if(d)printf("(min,max)=(%.2lf,%.2lf)\n",min,max);
    if(max==min)return;
    for(y=0;y<mat->h;y++)
    {
        h=y*mat->w;
        for(x=0;x<mat->w;x++)
        {
            col=mat->data[h+x];
            map=(col-min)/(max-min);
            if(map<0||map>1)printf("WTF\n");
            c=(unsigned char)(map*255);
            gfx_color(c,c,c);
            for(nx=0;nx<4;nx++)
            {
                for(ny=0;ny<4;ny++)
                {
                    gfx_point(4*x+nx,4*y+ny);
                }
            }
        }
    }
}

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
            if(mat->data[h+x]*mat->data[h+x]>ZERO)return 0;
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
    printf("\tPerforming QR decomposition...\n");
    qr_decomp* QR=qr_decompose(mat);
    matrix* RQ=multiply(QR->R,QR->Q);
    matrix* Q=Matrix(QR->Q->w,QR->Q->h,QR->Q->data),*tmp;
    unsigned long int i=0;
    printf("\tIterating until it's triangular...\n");
    while(!istriangular(RQ))
    {
        fflush(stdout);
        freeMatrix(&QR->R);
        freeMatrix(&QR->Q);
        free(QR);
        showmat(RQ,0);
        printf("\t\t(Internal QR decomposition)(%lu)\r",i++);
        fflush(stdout);
        QR=qr_decompose(RQ);
        freeMatrix(&RQ);
        tmp=multiply(Q,QR->Q);
        freeMatrix(&Q);
        Q=tmp;
        RQ=multiply(QR->R,QR->Q);
    }
    printf("\nFinished QR eigenpair algorithm.\n");
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
    ret=malloc(sizeof(SVD));
    unsigned long int non_null,i;
    printf("Generating square matrix...\n");
    matrix *trans=transposed(mat),*sq=multiply(trans,mat);
    showmat(sq,0);
    gfx_wait();
    printf("Calculating eigenpairs...\n");
    eigenpairs* eigen=qr_eigenpairs_square(sq);
    non_null=0;
    for(non_null=0;non_null<sq->w;non_null++)if(eigen->eigenvalues[non_null]<=ZERO)break;
    ret->VT=transposed_cut(eigen->eigenvectors,non_null);
    ret->U=Matrix(non_null,mat->h,NULL);
    /*For each eigenvalue*/
    for(i=0;i<non_null;i++)
    {
        eigen->eigenvalues[i]=sqrt(eigen->eigenvalues[i]);
        if(eigen->eigenvalues[i]>ZERO)
        {
            multiplyByVectToVect(mat,eigen->eigenvectors,ret->U,i,i);
            scaleCol(&ret->U,i,1/eigen->eigenvalues[i]);
        }
    }
    ret->S=square_from_vect(eigen->eigenvalues,non_null);
    ret->len=non_null;
    return ret;
}

matrix* reconstruct_from_n_values(SVD* svd,unsigned int n)
{
    if(!svd||n==0)return NULL;
    matrix* ret=Matrix(svd->VT->w,svd->U->h,NULL);
    unsigned int i,j,j2;
    unsigned long int hj,hj2,hr,eigen;
    for(i=0;i<n;i++)
    {
        eigen=svd->S->w*i+i;
        for(j=0;j<svd->U->h;j++)
        {
            hj=svd->U->w*j;
            hr=svd->VT->w*j;
            for(j2=0;j2<svd->VT->w;j2++)
            {
                hj2=svd->VT->w*i;
                ret->data[hr+j2]+=svd->S->data[eigen]*svd->U->data[hj+i]*svd->VT->data[hj2+j2];
            }
        }
    }
    return ret;
}