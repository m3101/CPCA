#include <stdlib.h>
#include <stdio.h>
#include "gfx.h"
#include "../src/cpca.h"
#include "../src/linal.h"

int main()
{
    matrix* mat=Matrix(10,10,NULL);
    int x,y,ly;
    for(y=0;y<mat->h;y++)
    {
        ly=mat->w*y;
        for(x=0;x<mat->w;x++)
            mat->data[ly+x]=x>=y?ly+x:0;
    }
    for(y=0;y<mat->h;y++)
    {
        ly=mat->w*y;
        for(x=0;x<mat->w;x++)
            printf("%.2lf\t",mat->data[ly+x]);
        printf("\n");
    }
    printf("###\n");
    sortMat(&mat);
    printmat(mat);
    printf("GAUSS:\n");
    row_gaussElim(&mat);
    printf("\nResult:\n");
    printmat(mat);
    return 0;
}