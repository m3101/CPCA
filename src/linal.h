#ifndef LINAL_H
#define LINAL_H
#define ZERO 0.00000000000001
/*
A simple collection of linear algebra functions
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

/*A matrix*/
typedef struct matrix{
    unsigned int w;/*Number of columns*/
    unsigned int h;/*Number of lines*/
    double* data;/*entries on the matrix (concatenated lines)*/
}matrix;

/*
Allocates a new matrix. If data!=NULL, also copies data from <data> to matrix.data
*/
matrix* Matrix(unsigned long int w, unsigned long int h, double* data);

/*
Allocates a new matrix using data from a char array (takes only 1 channel).
*/
matrix* Matrix_From_Img(unsigned long int w, unsigned long int h, char* data,char channel);

/*
Frees a matrix's allocated memory
*/
void freeMatrix(matrix** m);

/*
Scales a matrix's <n>th line by a factor of <s>
*/
void scaleLine(matrix** mat,unsigned long int n,double s);

/*
Scales a matrix's <n>th column by a factor of <s>
*/
void scaleCol(matrix** mat,unsigned long int n,double s);

/*
Swaps two of a matrix's lines for each other
*/
void swapLine(matrix** mat,unsigned long int n,unsigned long int m);

/*
Swaps two of a matrix's columns for one another
*/
void swapCol(matrix** mat,unsigned long int n,unsigned long int m);

/*
Sums a matrix's <m>th line's entries to its <n>th's
*/
void sumLines(matrix** mat,unsigned long int n, unsigned long int m);

/*
Sums a matrix's <m>th column's entries times <s> to its <n>th's
*/
void sumSCols(matrix** mat,unsigned long int n, unsigned long int m,double s);

/*
Sums a matrix's <m>th line's entries times <s> to its <n>th's
*/
void sumSLines(matrix** mat,unsigned long int n, unsigned long int m,double s);

/*
Transposes a matrix across its main diagonal
*/
matrix* transposed(matrix* mat);

/*
Returns a transposed matrix (across the main diagonal) ignoring elements after the <n>th column
*/
matrix* transposed_cut(matrix* mat,unsigned long int n);

/*
Multiplies two matrices
*/
matrix* multiply(matrix* a,matrix* b);

/*
Multiplies matrix <a> by <b>'s <v>th column
*/
matrix* multiplyByVect(matrix* a,matrix* b,unsigned long int v);

/*
Multiplies matrix <a> by <b>'s <v>th column and places the result in <dest>'s <vi>th column
*/
void multiplyByVectToVect(matrix* a,matrix* b,matrix* dest,unsigned long int v,unsigned long int vi);

/*
Sorts a matrix's lines so that the rightmost pivots are at the bottom
AND same-column pivots are ordered in descending order (top-to-bottom)
*/
void sortMat(matrix** mat);

/*
Brings a matrix into row echelon form
*/
void row_gaussElim(matrix ** mat);

/*
Prints a matrix
*/
void printmat(matrix* mat);

/*
A structure for storing QR-decompositions
*/
typedef struct qr_decomp{
    matrix* Q;
    matrix* R;
}qr_decomp;

/*
Applies the gram-schmidt process to the columns of a matrix
*/
void gram_schmidt_columns(matrix* mat);

/*
Performs qr-decomposition on a matrix
*/
qr_decomp* qr_decompose(matrix* mat);
#endif