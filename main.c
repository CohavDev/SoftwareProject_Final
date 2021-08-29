#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double **buildMatrixW(double **, int, int);
double **computeMatrixD(double **,int);
double calcWeight(double*, double*, int);
double** buildMatrixT(double **, int, int);
double * buildMatrixP(double **, int);
void diagonalStep(double **, double *, int);
void rotationMultiply(double **, double*,int);
double **buildPfromArray(double *, int);
double ** jacobiAlgorithm(double **, int);
void printMatrix(double**, int,int);

void mainFunction(int argc, char *argv){
    //TODO:continue method
    int k = (int)argv[1];
    if(k == 0){
        //TODO:call eigengap heuristic
    }

    char *goal = argv[2];
    if(strcmp(goal, "spk") == 0){
        //spectralKmeans();
    }
    if(strcmp(goal, "wam") == 0){
//        double **matW = buildMatrixW()
//        printMatrix(matW,n,n)
    }
    if(strcmp(goal, "ddg") == 0){
//        double **matD = computeMatrixD(points, n);
//            printMatrix(matD, n, n);
    }
    if(strcmp(goal, "lnorm") == 0){
        //print lnorm
    }
    if(strcmp(goal, "jacobi") == 0){
        //jacobiAlgorithm()
    }

}
void spectralkmeans(double** points, int n, int d){
    double **matW = buildMatrixW(points, d, n);
    //TODO:continue here
}
int main() {
    //testing
    int i,j;
    double arr[10][2] = {
            {-5.056, 11.011},
            {-6.409, -7.962},
            {5.694, 9.606},
            {6.606, 9.396},
            {-6.772, -5.727},
            {-4.498, 8.399},
            {-4.985, 9.076},
            {4.424, 8.819},
            {-7.595, -7.211},
            {-4.198, 8.371}};
    double **points = malloc(10*sizeof(double));
    for(i=0;i<10;i++){
        points[i] = malloc(2*sizeof (double));
        for(j=0;j<2;j++){
            points[i][j] = arr[i][j];
        }
    }
    //the weighted matrix test
    double **matW = buildMatrixW(points,2,10);
    printf("\nMatrix W:\n");
    printMatrix(matW,10,10);
    //the diagonal degree matrix test
    double **matD = computeMatrixD(matW,10);
    printf("\nMatrix D:\n");
    printMatrix(matD, 10,10);

    return 1;
}
//this function builds the weighted matrix W
double **buildMatrixW(double **points, int d,int n){
    int i, j;
    double tmp;
    //allocate memory for matrix Wnxn
    double** w = malloc(sizeof(double*)*n);
    for(i=0;i<n;i++){
        w[i] = malloc(n*sizeof(double));
    }
    //calculate matrix W
    for(i=0;i<n;i++){
        //each loop calculates row of upper diagonal and copies it to bottom
        w[i][i] = 0;
        for(j=i+1;j<n;j++){
            tmp = calcWeight(points[i],points[j], d);//TODO:check calcweight() correctness
            w[i][j] = tmp;
            w[j][i] = tmp; // matrix is symmetric
            //TODO: anyway for complexity improvements?
        }
    }
    return w;
}
//this function compute the diagonal matrix D
double **computeMatrixD(double **w, int n){
    //TODO: anyway for complexity improvements?
    int i, j;
    double tmp;
    //allocate memory for matrix Wnxn
    double** d = malloc(sizeof(double*)*n);
    for(i=0;i<n;i++){
        d[i] = calloc(n,sizeof(double));
    }
    for(i=0;i<n;i++){
        tmp = 0;
        for(j=0;j<n;j++){
            tmp+=w[i][j];
        }
        d[i][i] = tmp;
    }
    return d;
}

double calcWeight(double *a, double*b, int coordnum){
    //l2 norm calculation
    double sum = 0;
    int i;
    for(i=0;i<coordnum;i++){
        sum += (a[i]-b[i])*(a[i]-b[i]);
    }
    sum = sqrt(sum);
    return exp(-0.5*sum);

}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    int L[n1], R[n2];

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(int arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}
//build the normalized matrix T
double** buildMatrixT(double **matrixU, int n, int k){
    //allocate memory for matrix Tnxk
    int i,j;
    double sum;
    double** t = malloc(sizeof(double*)*n);
    for(i=0;i<n;i++){
        t[i] = malloc(k*sizeof(double));
    }
    //calc the 'normalizer' for each row
    double *norm = malloc(n*sizeof(double));
    for(i=0;i<n;i++){
        sum = 0;
        for(j=0;j<k;j++){
            sum += pow(matrixU[i][j],2);
        }
        norm[i] = sqrt(sum);
    }
    //calc matrix T
    for(i=0;i<n;i++){
        for(j=0;j<k;j++){
            t[i][j] = matrixU[i][j] / norm[i];
        }
    }
    return t;
}

//the method returns an array which contains: i,j,c,s. in this order
// a is a symmetric matrix of size nxn
//which represent matrix P
double * buildMatrixP(double **a, int n){
    //find Aij the maximum off-diagonal element of matrix A
    int l,m,i,j;
    double max = fabs(a[0][0]);
    i=0;
    j=0;
    for(l=0;l<n;l++){
        for(m=l+1;m++;m<n)
            if(fabs(a[l][m])>max){
                i=l;
                j=m;
                max = fabs(a[i][j]);
            }
        }

    //calculate c,s,t
    double teta = (a[j][j]-a[i][i])/(2*a[i][j]);
    int sign = teta>=0 ? 1 : -1;
    double t = sign / (fabs(teta)+sqrt((teta*teta) + 1));
    double c = 1 / sqrt((t*t) + 1);
    double s = t*c;
    //arr represents matrix P
    double* arr = malloc(sizeof(double)*4);
    arr[0] = i;
    arr[1] = j;
    arr[2] = c;
    arr[3] = s;
    return arr;

}
//this method calculates A' = P^t*A*P
void diagonalStep(double **a, double *p, int n){
    int i,j,m;
    double c,s;
    double *rowi,*rowj;//row i and row j of A'
    //for convenient
    i = (int)p[0];
    j = (int)p[1];
    c =p[2];
    s=p[3];
    //memory allocation
    rowi = calloc(n,sizeof(double));
    rowj = calloc(n,sizeof(double));

    //calc rowi and rowj
    for(m=0;m<n;m++){
        if(m == i){
            rowi[i] = (pow(c,2)*a[i][i]) + (pow(s,2)*a[j][j]) - (2*s*c*a[i][j]);
            rowj[i] = 0;
        }
        else if(m==j){
            rowj[m] = (pow(s,2)*a[i][i]) + (pow(c,2)*a[j][j]) - (2*s*c*a[i][j]);
            rowi[j] = 0;
        }
        else{
            rowi[m] = (c*a[m][i]) -(s*a[m][j]);
            rowj[m] = (c*a[m][j]) + (s*a[m][i]);
        }

    }
    //update matrix A to A'
    //update row i-th and j-th
    a[i] = rowi;
    a[j] = rowj;
    //update columns i-th and j-th. Remember A is symmetric
    for(m=0;m<n;m++){
        a[m][j] = rowj[m];
        a[m][i] = rowi[m];
    }
}
//calculates A*P
//the only changes in A are i-th and j-th columns
//TODO: any optimization possible? maybe A transpose?
void rotationMultiply(double **a, double* p,int n){
    int m,l,i,j;
    double coli,colj;
    i = (int)p[0];
    j=(int)p[1];
    for(m=0;m<n;m++){
        //p[2] = c
        //p[3] = s
        coli = a[m][i]*p[2] - a[m][j]*p[3];
        colj = a[m][i]*p[3] + a[m][j]*p[2];
        a[m][i] = coli;
        a[m][j] = colj;
    }
}
//method build 2-d matrix P from array that represent P
double **buildPfromArray(double *p, int n){
    int i,j,coli,colj;
    double **mat = malloc(sizeof(double*)*n);
    coli = (int)p[0];
    colj = (int)p[1];
    for(i=0;i<n;i++){
        mat[i] = calloc(n,sizeof(double ));
        mat[i][i] = 1;
        if(i == coli){
            mat[i][i] = p[2];//c
            mat[i][colj] = p[3];//s
        }
        else if(i == colj){
            mat[i][coli] = -p[3];//-s
            mat[i][colj] = p[2];//c
        }
    }
    return mat;
}
double ** jacobiAlgorithm(double **a, int n){
    //TODO: not done
    int converged = 0;
    double *p = buildMatrixP(a, n);//p1
    double **matV = buildPfromArray(p,n);//matrix V init (V = p1)
    diagonalStep(a, p, n);//calc A = A'
    //TODO:check convergence here too
    while(!converged){
        p = buildMatrixP(a, n); //pi
        rotationMultiply(matV, p, n);//V*pi
        diagonalStep(a,p,n);//calc A'
        converged = 1;//TODO:need to update according to convergence rule
    }
    return NULL;//TODO: need to output/return eigenvalues and eigenvectors
}
void printMatrix(double **a, int rows,int columns){
    int i,j;
    for(i=0;i<rows;i++){
        if(i!=0){
            printf("\n");
        }
        for(j=0;j<columns-1;j++){
            printf("%.3lf,",a[i][j]);//TODO:change to .4f?
        }
        printf("%.3lf",a[i][columns-1]);
    }
}
//double multiply(double**a, double **b,int r1,int c1, int r2,int c2){
//    assert(r1 == c2); // number of a' rows equals number of  b's columns
//    int i,j;
//    for(i=0;i<)
//}
