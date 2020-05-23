#include <stdlib.h>
#include <stdio.h>
#include "gfx.h"
#include "../src/cpca.h"
#include "../src/linal.h"

int main()
{
    double mm[]=
    {
        1,1,1,0,0,
        3,3,3,0,0,
        4,4,4,0,0,
        5,5,5,0,0,
        0,2,0,4,4,
        0,0,0,5,5,
        0,1,0,2,2
    };
    matrix* mat=Matrix(5,7,mm),*tmp,*prod;
    
    /*double mm[]=
    {
        1,1,
        0,1,
        1,0
    };
    matrix* mat=Matrix(2,3,mm),*tmp,*prod;*/
    printf("###\n");
    printmat(mat);
    SVD* svd=qr_svd(mat);
    tmp=multiply(svd->U,svd->S);
    prod=multiply(tmp,svd->VT);
    freeMatrix(&tmp);
    printmat(prod);
    return 0;
}