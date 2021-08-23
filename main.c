#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double **buildMatrixW(double **);
double **computeMatrixD(double **);
double calcWeight(double*, double*);
int n ,k,coordnum;
int main() {
    printf("Hello, World!\n");
    return 0;
}
//this function builds the weighted matrix W
double **buildMatrixW(double **points){
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
        for(j=i;j<n;j++){
            tmp = calcWeight(points[i],points[j]);//TODO:check calcweight() correctness
            w[i][j] = tmp;
            w[j][i] = tmp; // matrix is symmetric
            //TODO: anyway for complexity improvements?
        }
    }
    return w;
}
//this function compute the diagonal matrix D
double **computeMatrixD(double **w){
    //TODO:should create method for matrix allocation?
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

double calcWeight(double *a, double*b){
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
    //calc the normalizer for each row
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
//double multiply(double**a, double **b,int r1,int c1, int r2,int c2){
//    assert(r1 == c2); // number of a' rows equals number of  b's columns
//    int i,j;
//    for(i=0;i<)
//}
