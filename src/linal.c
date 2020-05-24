#define _GNU_SOURCE
#include "linal.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*
A simple collection of naïve linear algebra functions
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
Allocates a new matrix. If data!=NULL, also copies data from <data> to matrix.data
*/
matrix* Matrix(unsigned long int w, unsigned long int h, double* data)
{
    if(w==0||h==0)return NULL;
    matrix* ret=malloc(sizeof(matrix));
    unsigned long int len=w*h;
    ret->w=w;
    ret->h=h;
    ret->data=calloc(w*h,sizeof(double));
    if(data)
    {
        for(unsigned long int i=0;i<len;i++) ret->data[i]=data[i];
    }
    return ret;
}

/*
Allocates a new matrix using data from a char array (takes only 1 channel).
*/
matrix* Matrix_From_Img(unsigned long int w, unsigned long int h, char* data,char channel)
{
    if(w==0||h==0)return NULL;
    matrix* ret=malloc(sizeof(matrix));
    unsigned long int len=w*h;
    ret->w=w;
    ret->h=h;
    ret->data=calloc(w*h,sizeof(double));
    if(data)
    {
        for(unsigned long int i=0;i<len;i++) ret->data[i]=(unsigned char)data[(i*4)+channel];
    }
    return ret;
}

/*
Frees a matrix's allocated memory
*/
void freeMatrix(matrix** m)
{
    if(!m||!*m)return;
    free((*m)->data);
    free(*m);
    *m=NULL;
}

/*
Scales a matrix's <n>th line by a factor of <s>
*/
void scaleLine(matrix** mat,unsigned long int n,double s)
{
    unsigned long int h;
    if(!mat||!*mat||n>=(*mat)->h)return;
    h=(*mat)->w*n;
    for(unsigned long int i=0;i<(*mat)->w;i++)(*mat)->data[i+h]*=s;
}
void scaleCol(matrix** mat,unsigned long int n,double s)
{
    if(!mat||!*mat||n>=(*mat)->w)return;
    for(unsigned long int i=0;i<(*mat)->h;i++)
        (*mat)->data[n+(i*(*mat)->w)]*=s;
}

/*
Swaps two of a matrix's lines for each other
*/
void swapLine(matrix** mat,unsigned long int n,unsigned long int m)
{
    double tmp;
    unsigned long int h1,h2;
    if(!mat||!*mat||m>=(*mat)->h||n>=(*mat)->h)return;
    h1=(*mat)->w*m;
    h2=(*mat)->w*n;
    for(unsigned long int i=0;i<(*mat)->w;i++)
    {
        tmp=(*mat)->data[i+h1];
        (*mat)->data[i+h1]=(*mat)->data[i+h2];
        (*mat)->data[i+h2]=tmp;
    }
}

/*
Swaps two of a matrix's columns for one another
*/
void swapCol(matrix** mat,unsigned long int n,unsigned long int m)
{
    double tmp;
    unsigned long int h1;
    if(!mat||!*mat||m>=(*mat)->w||n>=(*mat)->w)return;
    for(unsigned long int i=0;i<(*mat)->h;i++)
    {
        h1=(*mat)->w*i;
        tmp=(*mat)->data[h1+n];
        (*mat)->data[h1+n]=(*mat)->data[h1+m];
        (*mat)->data[h1+m]=tmp;
    }
}

/*
Sums a matrix's <m>th line's entries to its <n>th's
*/
void sumLines(matrix** mat,unsigned long int n, unsigned long int m)
{
    unsigned long int h1,h2;
    if(!mat||!*mat||m>=(*mat)->h||n>=(*mat)->h)return;
    h1=(*mat)->w*m;
    h2=(*mat)->w*n;
    for(unsigned long int i=0;i<(*mat)->w;i++)
        (*mat)->data[i+h2]+=(*mat)->data[i+h1];
}

/*
Sums a matrix's <m>th line's entries times <s> to its <n>th's
*/
void sumSLines(matrix** mat,unsigned long int n, unsigned long int m,double s)
{
    unsigned long int h1,h2;
    if(s==0||!mat||!*mat||m>=(*mat)->h||n>=(*mat)->h)return;
    h1=(*mat)->w*m;
    h2=(*mat)->w*n;
    for(unsigned long int i=0;i<(*mat)->w;i++)
        (*mat)->data[i+h2]+=(*mat)->data[i+h1]*s;
}

void sumSCols(matrix** mat,unsigned long int n, unsigned long int m,double s)
{
    unsigned long int h1;
    if(s==0||!mat||!*mat||m>=(*mat)->w||n>=(*mat)->w)return;
    for(unsigned long int i=0;i<(*mat)->h;i++)
    {
        h1=i*(*mat)->w;
        (*mat)->data[n+h1]+=(*mat)->data[m+h1]*s;
    }
}
/*
Returns a transposed matrix (across the main diagonal)
*/
matrix* transposed(matrix* mat)
{
    if(!mat)return NULL;
    matrix* ret=Matrix(mat->h,mat->w,NULL);
    unsigned long int x,y,ly;
    for(y=0;y<mat->h;y++)
    {
        ly=y*mat->w;
        for(x=0;x<mat->w;x++)
            ret->data[(x*mat->h)+y]=mat->data[ly+x];
    }
    return ret;
}

/*
Returns a transposed matrix (across the main diagonal) ignoring elements after the <n>th column
*/
matrix* transposed_cut(matrix* mat,unsigned long int n)
{
    if(!mat)return NULL;
    matrix* ret=Matrix(mat->h,n,NULL);
    unsigned long int x,y,ly;
    for(y=0;y<mat->h;y++)
    {
        ly=y*mat->w;
        for(x=0;x<n;x++)
            ret->data[(x*mat->h)+y]=mat->data[ly+x];
    }
    return ret;
}

/*Returns the index of the first non-zero element on a line*/
unsigned long int pivot(double* line,unsigned int w)
{
    unsigned long int x;
    for(x=0;x<w&&line[x]==0;x++);
    return x;
}

/*
Multiplies two matrices
*/
matrix* multiply(matrix* a,matrix* b)
{
    if(a->w!=b->h)return NULL;
    matrix* ret=Matrix(b->w,a->h,NULL);
    unsigned long int x,y,i,lya,lyb,lyr;
    for(y=0;y<a->h;y++)
    {
        lya=y*a->w;
        lyr=y*b->w;
        for(x=0;x<b->w;x++)
        {
            for(i=0;i<a->w;i++)
            {
                lyb=i*b->w;
                ret->data[lyr+x]+=a->data[lya+i]*b->data[lyb+x];
            }
        }
    }
    return ret;
}

/*
Multiplies matrix <a> by <b>'s <v>th column
*/
matrix* multiplyByVect(matrix* a,matrix* b,unsigned long int v)
{
    if(a->w!=b->h)return NULL;
    matrix* ret=Matrix(1,a->h,NULL);
    unsigned long int y,i,lya,lyb;
    for(y=0;y<a->h;y++)
    {
        lya=y*a->w;
        for(i=0;i<a->w;i++)
        {
            lyb=i*b->w;
            ret->data[y]+=a->data[lya+i]*b->data[lyb+v];
        }
    }
    return ret;
}

/*
Multiplies matrix <a> by <b>'s <v>th column and places the result in <dest>'s <vi>th column
*/
void multiplyByVectToVect(matrix* a,matrix* b,matrix* dest,unsigned long int v,unsigned long int vi)
{
    if(a->w!=b->h||dest->h!=a->h||vi>=dest->w)return;
    unsigned long int y,i,lya,lyb,lyd;
    for(y=0;y<a->h;y++)
    {
        lya=y*a->w;
        lyd=y*dest->w;
        for(i=0;i<a->w;i++)
        {
            lyb=i*b->w;
            dest->data[lyd+vi]+=a->data[lya+i]*b->data[lyb+v];
        }
    }
}

int sortMat_compar(const void* _a,const void* _b,void* arg)
{
    const double* a=_a,*b=_b;
    unsigned long int w=*(unsigned long int*)arg;
    unsigned long int pivotA,pivotB;
    for(pivotA=0;pivotA<w&&a[pivotA]==0;pivotA++);
    for(pivotB=0;pivotB<w&&b[pivotB]==0;pivotB++);
    if(pivotB==pivotA)
    {
        if(pivotA==w)return 0;
        return a[pivotA]>b[pivotB]?1:-1;
    }
    else
        return pivotA>pivotB?1:-1;
}

/*
Sorts a matrix's lines so that the rightmost pivots are at the bottom
AND same-column pivots are ordered in descending order (top-to-bottom)
*/
void sortMat(matrix** mat)
{
    if(!mat||!(*mat))return;
    qsort_r((*mat)->data,(*mat)->h,(*mat)->w*sizeof(double),sortMat_compar,&(*mat)->w);
}

/*
Brings a matrix into row echelon form
*/
void row_gaussElim(matrix ** mat)
{
    if(!mat||!(*mat))return;
    sortMat(mat);
    unsigned long int x,y,ly,y2,ly2;
    for(y=0;y<(*mat)->h;y++)
    {
        ly=(*mat)->w*y;
        x=pivot(&(*mat)->data[ly],(*mat)->w);
        if(x==(*mat)->w)return;
        for(y2=y+1;y2<(*mat)->h;y2++)
        {
            ly2=(*mat)->w*y2;
            sumSLines(mat,y2,y,-(*mat)->data[ly2+x]/(*mat)->data[ly+x]);
        }
        scaleLine(mat,y,1/(*mat)->data[ly+x]);
    }
}

void printmat(matrix* mat)
{
    unsigned long int x,y,ly;
    for(y=0;y<mat->h;y++)
    {
        ly=mat->w*y;
        for(x=0;x<mat->w;x++)
            printf("%.4lf\t",mat->data[ly+x]);
        printf("\n");
    }
}

double scalProd(matrix* mat,unsigned long int c1,unsigned long int c2)
{
    unsigned long int y,jump;
    double ret=0;
    if(!mat)return 0;
    for(y=0;y<mat->h;y++)
    {
        jump=y*mat->w;
        ret+=mat->data[jump+c1]*mat->data[jump+c2];
    }
    return ret;
}

void gram_schmidt_columns(matrix* mat)
{
    if(!mat)return;
    unsigned long int x,x1,r,i;
    double s;
    for(x=0;x<mat->w;x++)
    {
        i=0;
        s=scalProd(mat,x,x);
        while(s<=ZERO)
        {
            r=(mat->w-x)-i;
            if(r==0)return;
            for(x1=x;x1<mat->w-1;x1++)
                swapCol(&mat,x1,x1+1);
            i++;
            s=scalProd(mat,x,x);
        }
        for(x1=0;x1<x;x1++)
        {
            s=scalProd(mat,x1,x1);
            if(s==0) continue;
            else s=scalProd(mat,x,x1)/s;
            sumSCols(&mat,x,x1,-s);
        }
        s=sqrt(scalProd(mat,x,x));
        if(s<=ZERO&&x<mat->w-1)
        {
            for(x1=x;x1<mat->w-1;x1++)
                swapCol(&mat,x1,x1+1);
            x--;
        }
        else
            scaleCol(&mat,x,1/s);
    }
}

qr_decomp* qr_decompose(matrix* mat)
{
    if(!mat)return NULL;
    qr_decomp* ret=malloc(sizeof(qr_decomp));
    matrix* trans;
    ret->Q=Matrix(mat->w,mat->h,mat->data);
    gram_schmidt_columns(ret->Q);
    trans=transposed(ret->Q);
    ret->R=multiply(trans,mat);
    freeMatrix(&trans);
    return ret;
}