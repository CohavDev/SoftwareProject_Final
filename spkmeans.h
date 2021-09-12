#ifndef STRUCTEIGEN
#define STRUCTEIGEN
typedef struct Eigen{
    int indexOfVector;
    double eigenValue;
}Eigen;
typedef struct Cluster {
    int count;
    /**pointers to arrays**/
    double *sum;
    double *mean;

}Cluster;
#endif
double **buildMatrixW(double **, int, int);
double **computeMatrixD(double **,int);
double calcWeight(double*, double*, int);
void multDiagonalLeft(double**, double**, int);
void multDiagonalRight(double**, double**, int);
double **buildMatrixLnorm(double **points, int d, int n);
double** buildMatrixT(double **, int, int);
void buildMatrixP(double **, int,double*);
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
void freeMem(double** vecArray, Cluster** clusterArray,int n,int k);
void kMeans(double** vecArray,Cluster** clusterArray, int n, int k);
void printJacobi(double** matrix, int n);
Cluster ** initClusters(double **vecArray,int n, int k);
/**kmeans - HW1 **/
int reCalcMeans();
int calcMean(Cluster *clust,int d);
void addVector(Cluster *clust, double *vec, int d);
void findCluster(double *vec, struct Cluster** clusterArray,int k);
double distance(double *x, double *y,int d);
void refreshClusters(Cluster** clusterArray,int k,int d);
void printMeans(Cluster **clusterArray,int k);
void initFromFile(int k, char* fileName, char* goal);