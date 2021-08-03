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
//double multiply(double**a, double **b,int r1,int c1, int r2,int c2){
//    assert(r1 == c2); // number of a' rows equals number of  b's columns
//    int i,j;
//    for(i=0;i<)
//}
