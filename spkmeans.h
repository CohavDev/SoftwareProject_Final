#ifndef STRUCTEIGEN
#define STRUCTEIGEN
typedef struct Eigen{
    int indexOfVector;
    double eigenValue;
}Eigen;
#endif
double **buildMatrixW(double **, int, int);
double **computeMatrixD(double **,int);
double calcWeight(double*, double*, int);
void multDiagonalLeft(double**, double**, int);
void multDiagonalRight(double**, double**, int);
double **buildMatrixLnorm(double **points, int d, int n);
double** buildMatrixT(double **, int, int);
double * buildMatrixP(double **, int);
double diagonalStep(double **, double *, int, double);
void rotationMultiply(double **, double*,int);
double **buildPfromArray(double *, int);
double calcConvergenceValue(double**, int);
double ** jacobiAlgorithm(double **, int);
void printMatrix(double**, int,int);
double** handleGoal(double **, int , int , int, char *, int);
int compare(const void *a, const void *b);
double** buildMatrixU(double** matV, int n, int k);
int determineK(struct Eigen* structArray, int n);