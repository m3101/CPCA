#include <stdlib.h>
#include <stdio.h>
#include "gfx.h"
#include "../src/cpca.h"
#include "../src/linal.h"
#include "../lib/bmpreader.h"

int main()
{
    int w,h,i,prec=0,j;
    double eigenweight,accweight;
    char* img=readbmp("satie.bmp"),c;
    hwfromsize(sizebmp("satie.bmp"),w,h);
    matrix* mat=Matrix_From_Img(w,h,img,1);
    gfx_open(w*4,h*4,"SVD!");
    showmat(mat,0);
    gfx_wait();
    printf("Starting SVD...\n");
    SVD* svd=qr_svd(mat);
    printf("Done!\nReconstructing...\n");
    matrix* reconst=reconstruct_from_n_values(svd,svd->len-prec);
    printf("Eigenvalues: {");
    eigenweight=0;
    for(i=0;i<svd->len;i++)
    {printf("%.2lf,",svd->S->data[i+i*svd->S->w]);eigenweight+=svd->S->data[i+i*svd->S->w];}
    printf("0}\n");
    showmat(reconst,0);
    gfx_wait();
    while((c=gfx_wait())!='q')
    {
        freeMatrix(&reconst);
        switch (c)
        {
        case 'w':
            if(prec>1)
                prec--;
            break;
        case 's':
            if(prec<svd->len-1)
                prec++;
            break;
        default:
            break;
        }
        accweight=0;
        for(j=0;j<svd->len-prec;j++)accweight+=svd->S->data[j+j*svd->S->w];
        reconst=reconstruct_from_n_values(svd,svd->len-prec);
        printf("Accumulated 'eigenweight': %.2lf\r",accweight/eigenweight);
        fflush(stdout);
        showmat(reconst,0);
    }
    return 0;
}