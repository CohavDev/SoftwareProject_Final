#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "spkmeans.h"

void spectralKmeans(double **points, int n, int d, int k){
    double** mat;
    Cluster **clusterArray;
    int i;
    mat = buildMatrixLnorm(points,d,n);
    mat = jacobiAlgorithm(mat, n);
    mat = buildMatrixU(mat,n,k);
    //compute matrix T
    mat = buildMatrixT(mat, n,k);
    //call original Kmeans algorithm from HW1 on mat(matrix T) TODO:finish the call and print
    clusterArray = initClusters(mat, n, k);
    kMeans(mat,clusterArray,n,k);
    printMeans(clusterArray,k);
    freeMem(mat,clusterArray,n,k);//free mat and clusterArray


}
void initCmd(char* argv[]){
    int k;
    char *goal, *fileName;
    k = strtol(argv[1],NULL,10);
    goal = argv[2];
    fileName = argv[3];
    initFromFile(k,fileName,goal);
}
void initFromFile(int k, char* fileName,char *goal){
    /*find paramater d with first line*/
    int n,d;
    char *str;
    double number;
    int i,j,size,reached_end;
    double *arr;
    double **vecArray;
    const int LINE_MAX_LENGTH = 1000; /* according to forum*/
    char* line;
    char buffer[LINE_MAX_LENGTH];
    FILE *fp = fopen(fileName,"r");
    if(fp == NULL){
        printf("\nAn Error Has Occured\n");
        assert(fp!=NULL);
    }
    d = 1;
    n = 0;
    //gets d , num of coordinates of each vector/point
    line = fgets(buffer,LINE_MAX_LENGTH,fp);
    str = strtok(line,",");
    while(str!=NULL){
        str = strtok(NULL,",");
        d++;
    }
    rewind(fp);
    size = 50;
    vecArray = malloc(sizeof(double*) * size); /*init to size of size*/
    assert(vecArray !=NULL);
    reached_end = 0; /*"False"*/

    /*Read vectors from file*/
    while (reached_end == 0) {
        arr  = malloc(sizeof(double) * d);  /*array of coordinates*/
        if(n == size){
            vecArray = realloc(vecArray, (2*size*sizeof(double*)));
            assert(vecArray !=NULL);
            size += size;
        }
        line = fgets(buffer,LINE_MAX_LENGTH,fp);
        str = strtok(line,",");
        for (i = 0; i < d; i++) {
            number = atof(str);
            if(feof(fp)){
                reached_end = 1;
                break;
            }
            arr[i] = number;
            str = strtok(NULL, ",");

        }
        if(reached_end != 1){
            n++;
            vecArray[n-1] = arr;
        }

    }
    assert(k >= 0 && k <= n);
    vecArray = realloc(vecArray, sizeof(double*) * n);
    assert(vecArray != NULL);
    handleGoal(vecArray,n,d,k,goal,1);
}
//sort the eigenvalues and its matching eigenvectors
//matV is (n+1) x n
double** buildMatrixU(double** matV, int n, int k){
    struct Eigen *structArray;
    int i,j, index;
    //init structArray
    structArray = malloc(n * sizeof(Eigen));
    for(i=0;i<n;i++){
        structArray[i].eigenValue = matV[n][i];
        structArray[i].indexOfVector = i;
    }
    //sort structArray with eigenvalue of each struct as a key
    qsort(structArray,n,sizeof(Eigen),compare);
//    printf("\nSorted EigenValues: \n");
//    for(i=0;i<n;i++){
//        printf("%f, ",structArray[i].eigenValue);
//    }
//    printf("\n");
    //determine k
    if(k == 0){
        k = determineK(structArray, n);
    }
    //build matrix U
    double **matrixU = malloc(n*sizeof(double*));

//    for(i=0;i<n;i++){
//        matrixU[i] = malloc(k*sizeof(double));
//    }
//    for(j=0;j<k;j++){
//        if(j==1){
//            printf("hey");
//        }
//        index = structArray[j].indexOfVector;
//        for(i=0;i<n;i++){
//            matrixU[i][j] = matV[i][index];
//        }
//    }
    for(i=0;i<n;i++){
        matrixU[i] = malloc(k*sizeof(double));
        for(j=0;j<k;j++){
            index = structArray[j].indexOfVector;
            matrixU[i][j] = matV[i][index];
        }
    }
    //free matrixV
    for(i=0;i<n+1;i++){
        free(matV[i]);
    }
    free(matV);
    //free sturctArray
    free(structArray);
    return matrixU;

}
//compare two Eigen Struct
int compare(const void *a, const void *b){
    Eigen *s1 = (Eigen*)a;
    Eigen *s2 = (Eigen*)b;
    double result = s1->eigenValue - s2->eigenValue;
    if(result == 0){
        return 0;
    }
    return (result>0) ? 1 : -1;
}
//TODO: should replace this this main()
//print: 1 - called from C. 0 - called from Python via module
double** handleGoal(double **points, int n, int d, int k, char *goal, int print){
    //TODO:continue method
    double **mat;
    int i;

    if(strcmp(goal, "spk") == 0){
        spectralKmeans(points,n,d,k);
        return NULL;//nothing to print or return here TODO:yes?
    }
    if(strcmp(goal, "wam") == 0){
        mat = buildMatrixW(points, d, n);
    }
    if(strcmp(goal, "ddg") == 0){
        mat = computeMatrixD(buildMatrixW(points,d,n), n);
    }
    if(strcmp(goal, "lnorm") == 0){
        mat = buildMatrixLnorm(points, d, n);
    }
    if(strcmp(goal, "jacobi") == 0){
        mat = jacobiAlgorithm(points,n);
    }
    if(print == 1){
        //called from C
        //free points
        for(i=0;i<n;i++){
            free(points[i]);
        }
        free(points);
    }

    //print or return
    if(print == 1){
        if(strcmp(goal, "jacobi") == 0){
            printJacobi(mat,n);
            //free
            for(i=0;i<n+1;i++){
                free(mat[i]);
            }
        }
        else{
            printMatrix(mat,n,n);
            //free
            for(i=0;i<n;i++){
                free(mat[i]);
            }
        }
        free(mat);
        return NULL;
    }
    else{
        //return to module
        return mat;
    }
}
void printJacobi(double** matrix, int n){
    //print eigenvalues
    int i,j;
    for(i=0;i<n-1;i++){//TODO: change to .4f
        printf("%.4f,",matrix[n][i]);
    }
    printf("%.4f\n",matrix[n][n-1]);
    //print eigenvectors
    for(j=0;j<n;j++){
        for(i=0;i<n-1;i++){
            printf("%.4f,",matrix[i][j]);
        }
        printf("%.4f\n",matrix[n-1][j]);
    }

}
int main(int argc, char* argv[]){
    initCmd(argv);
}
int test() {
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

    //lnorm matrix test
    double **matLNorm = buildMatrixLnorm(points, 2, 10);
    printf("\nMatrix Lnorm: \n");
    printMatrix(matLNorm, 10, 10);

    //jacobi test
    double** matJacobi = jacobiAlgorithm(matLNorm, 10);
    double** matrixU = buildMatrixU(matJacobi, 10, 0);
    printf("\nMatrix U:\n ");
    printMatrix(matrixU, 10,2);
    printf("\n");
    printf("\nMatrix T:\n ");
    printMatrix(buildMatrixT(matrixU,10,2), 10,2);
    printf("\n");
    printf("\ntest handleGoal: \n");
    handleGoal(points,10,2,2,"jacobi",1);
    handleGoal(points,10,2,2,"spk",1);
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
            tmp = calcWeight(points[i],points[j], d);
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
//multiply D*A
void multDiagonalLeft(double **diagonal, double **a, int n){
    int i, j;
    double multiplier;
    for(i=0;i<n;i++){
        multiplier = diagonal[i][i];
        for(j=0;j<n;j++){
            a[i][j] = a[i][j] * multiplier;
        }
    }
}
//multiply A*D
void multDiagonalRight(double **a, double **diagonal,int n) {

    int i,j;
    double multiplier;
    for (j=0; j<n;j++) {
        multiplier = diagonal[j][j];
        for (i=0;i<n;i++) {
            a[i][j] = a[i][j] * multiplier;
        }
    }

}

//calculate matrix Lnorm
double **buildMatrixLnorm(double **points, int d, int n) {
    double **matrixW = buildMatrixW(points, d, n);
    double **matrixD = computeMatrixD(matrixW, n);
    int i, j;
    double temp;
    //calc D^-0.5
    for (i = 0; i < n; i++) {
        temp = 1 / sqrt(matrixD[i][i]);
        matrixD[i][i] = temp;
    }
    multDiagonalLeft(matrixD, matrixW, n);
    multDiagonalRight(matrixW, matrixD, n);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            temp = matrixW[i][j];
            if(i==j){
                matrixW[i][j] = 1 - temp;
            }
            else{
                matrixW[i][j] = 0 - temp;
            }
        }
    }
    //free memory of matrixD
    for(i=0;i<n;i++){
        free(matrixD[i]);
    }
    free(matrixD);
    //return lnrom
    return matrixW;//lnorm was insert into matrixW, free is not here
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
    //free norm
    free(norm);
    //free matrixU
    for(i=0;i<n;i++){
        free(matrixU[i]);
    }
    free(matrixU);
    return t;
}
//the eigengap hueristic
//eigenvalues are sorted in rising sequence
int determineK(Eigen* structArray, int n) {
    int i, m;
    int argi;
    double max, temp;
    Eigen *eigen;
    max = 0;
    m = n/2;
    for (i = 0; i<=m; i++) {
        temp = structArray[i+1].eigenValue- structArray[i].eigenValue;
        if(temp >= max){//TODO: >= or >?
            max = temp;
            argi = i;
        }
    }
    return argi;
}

//the method saves to 'arr' an array which contains: i,j,c,s. in this order
// a is a symmetric matrix of size nxn
//which represent matrix P
void buildMatrixP(double **a, int n, double *arr){
    //find Aij the maximum off-diagonal element of matrix A
    int l,m,i,j;
    double max = fabs(a[0][1]);
    i=0;
    j=1;
    for(l=0;l<n;l++){
        for(m=l+1;m<n;m++)
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
//    double* arr = malloc(sizeof(double)*4);
    arr[0] = i;
    arr[1] = j;
    arr[2] = c;
    arr[3] = s;

}
//this method calculates A' = P^t*A*P
double diagonalStep(double **a, double *p, int n, double convergenceValueA){
    int i, j, r;
    double c, s, newConvergenceValue;
    double addToSum, subtractFromSum;
    double *coli,*colj;//row i and row j of A'
    //for convenient
    i = (int)p[0];
    j = (int)p[1];
    c =p[2];
    s=p[3];
    //memory allocation
    coli = calloc(n, sizeof(double));
    colj = calloc(n, sizeof(double));
    addToSum = 0;
    subtractFromSum = 0;
    //calc coli and colj
    for(r=0; r < n; r++){
        if(r == i){
            coli[i] = (pow(c, 2) * a[i][i]) + (pow(s, 2) * a[j][j]) - (2 * s * c * a[i][j]);//A'ii
            colj[i] = 0;
            subtractFromSum += 2*pow(a[j][i], 2);
        }
        else{
            if(r == j){
                colj[j] = (pow(s, 2) * a[i][i]) + (pow(c, 2) * a[j][j]) + (2 * s * c * a[i][j]);//A'jj
                coli[j] = 0;
            }
            else{
                coli[r] = (c * a[r][i]) - (s * a[r][j]);//A'ri
                colj[r] = (c * a[r][j]) + (s * a[r][i]);//A'rj
                //calc off(A')^2
                addToSum += 2*(pow(coli[r], 2) + pow(colj[r], 2));//counted twice for columns also
                subtractFromSum += 2 * (pow(a[i][r], 2) + pow(a[j][r], 2));//subtract A's elements
            }
        }
    }
    newConvergenceValue = convergenceValueA + addToSum - subtractFromSum;
    //update matrix A to A'
    //update row i-th and j-th
    a[i] = coli;
    a[j] = colj;
    //update columns i-th and j-th. Remember A is symmetric
    for(r=0; r < n; r++){
        a[r][j] = colj[r];
        a[r][i] = coli[r];
    }
    return newConvergenceValue; // return off(A')^2
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
//calculate convergence off(A)^2 if matrix a
double calcConvergenceValue(double **a, int n){
    double squareSum;
    int i, j;
    squareSum = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i!=j){
                squareSum+= pow(a[i][j], 2);
            }
        }
    }
    return squareSum;
}
//return 1 if converged , 0 otherwise
int isConverged(double n1 , double n2, double eps){
    if(n1 - n2 <= eps){
        return 1;
    }
    return 0;
}
double ** jacobiAlgorithm(double **a, int n){
    int converged;//0 False, 1 True
    int index;
    int iters = 0;//number of iterations
    double epsilon = pow(10,-15);
    double *p = malloc(sizeof(double)*4);
    buildMatrixP(a, n, p);//p1
    double convergenceValueA = calcConvergenceValue(a, n);
    double **matV = buildPfromArray(p,n);//matrix V init (V = p1)
//    printf("\nIteration No. 1\ni = %f, j = %f\n",p[0],p[1]);
    double newConvergenceValue = diagonalStep(a, p, n, convergenceValueA);//calc A = A'
//    printMatrix(a,n,n);
    converged = isConverged(convergenceValueA, newConvergenceValue, epsilon);//checking convergence
    iters++;
    while((converged == 0) && iters < 100){
        convergenceValueA = newConvergenceValue;
        buildMatrixP(a, n,p); //pi
//        printf("\nIteration No. %d\ni = %f,  j = %f\n",iters+1,p[0],p[1]);
        rotationMultiply(matV, p, n);//V*pi
        newConvergenceValue = diagonalStep(a,p,n, convergenceValueA);//calc A'
//        printf("\nMatrix A': \n");
//        printMatrix(a,n,n);
        converged = isConverged(convergenceValueA, newConvergenceValue, epsilon);//checking convergence
        iters++;
    }
//    printf("\nEigenValues: \n");
//    for(index = 0;index<n;index++){
//        printf("%.3f, ",a[index][index]);
//    }
//    printf("\nMatrix V:\n");
//    printMatrix(matV,n,n);
    free(p);
    //attach to matV the eigenvalues as the last row
    matV = realloc(matV, sizeof(double*)*(n+1));
    matV[n] = malloc(n*sizeof(double));
    for(index = 0;index<n;index++){
        matV[n][index] = a[index][index];
    }
//    printf("\nFinal Matrix of jacobi:\n");
//    printMatrix(matV, n+1,n);
    return matV;
}
void printMatrix(double **a, int rows,int columns){
    int i,j;
    for(i=0;i<rows;i++){
        if(i!=0){
            printf("\n");
        }
        for(j=0;j<columns-1;j++){
            if(a[i][j]<0 && a[i][j] > -0.0005){
                printf("%lf,",0.0000);
            }
            else{
                printf("%.4lf,",a[i][j]);
            }
        }
        if(a[i][columns-1]<0 && a[i][columns-1] > -0.0005){
            printf("%lf",0.0000);
        }
        else{
            printf("%.4lf",a[i][columns-1]);
        }
    }
}

/**paste HW1 code below here **/
/** global variables **/
//double **vecArray;
//Cluster **clusterArray;
const int maxIter = 300;

/** functions **/

/*input: Array of pointers of clusters and its length
output: None, calcs means of all clusters*/
int reCalcMeans(Cluster **clusterArray,int k) {
    int changed = 0;  /*"False"*/
    int i;
    for (i = 0; i < k; i++) {
        if (calcMean(clusterArray[i],k) == 1) {
            changed = 1;
        }
    }
    return changed;
}

/*calcs mean of specific cluster
returns 1 if changed, 0 otherwise*/
int calcMean(Cluster *clust,int d) {
    int i;
    int changed = 0; /*"False"*/
    double *sum, *mean; /*pointers to sum and mean of clust*/
    double calcVal;
    sum = clust->sum;
    mean = clust->mean;
    for (i = 0; i < d; i++) {
        calcVal = (sum[i] / (*clust).count);
        if (mean[i] != calcVal) {
            changed = 1; /*"True"*/
        }
        mean[i] = calcVal;
    }
    return changed;
}

/*adds specific vector to specific cluster*/
void addVector(Cluster *clust, double *vec, int d) {
    int i;
    double *sum;
    sum = (*clust).sum;
    for (i = 0; i < d; i++) {
        sum[i] = sum[i] + vec[i];
    }
    (*clust).count++;
}

/*input: Array of clusters, it's length and a vector.
output: None. the function find the closest cluster to the vector
 and insert vector to cluster.*/
void findCluster(double *vec, Cluster** clusterArray, int k) {
    /*declarations*/
    int i;
    Cluster *minCluster;
    double minDistance, tempDistance;

    /*default minimum*/
    minCluster = clusterArray[0];
    minDistance = distance(vec, minCluster->mean,k);

    /*loop through all clusters, except for 1st*/
    for (i = 1; i < k; i++) {
        tempDistance = distance(vec, (clusterArray[i])->mean, k);
        if (tempDistance < minDistance) {
            minDistance = tempDistance;
            minCluster = clusterArray[i];
        }
    }
    addVector(minCluster, vec,k);
}
/*init clusters*/
Cluster ** initClusters(double **vecArray,int n, int k){
    Cluster** clusterArray;
    int i,j;
    clusterArray = malloc(k * sizeof(Cluster*));
    assert(clusterArray != NULL);
    for(i=0; i < k;i++){
        Cluster *cl = malloc(sizeof(Cluster));
        cl->mean = malloc(k*sizeof(double));
        cl->sum = calloc(k,sizeof(double ));
        assert(cl->sum !=NULL && cl->mean != NULL);
        /*deep copy*/
        for(j=0;j<k;j++){
            (cl->mean)[j] = vecArray[i][j];
        }
        cl->count = 0;
        clusterArray[i] = cl;
    }
    return clusterArray;
}
//print the final centroids
void printMeans(Cluster **clusterArray,int k) {
    int i, j;
    double *temp;
    /*loop through clusters*/
    for (i = 0; i < k; i++) {
        temp = (clusterArray[i])->mean;
        /*loop through mean array*/
        for (j = 0; j < k-1; j++) {
            printf("%0.4f,", temp[j]);
        }
        printf("%0.4f\n", temp[k-1]);
    }
}
/*kMeans function*/
void kMeans(double** vecArray,Cluster** clusterArray, int n, int k) {
    int i, iterCount, changed;
    changed = 1;
    iterCount = 0;
//    initFromFile(k, clusterArray, vecArray);
    while (iterCount < maxIter && changed == 1) {
        /*loop through vectors and insert to clusters*/
        for (i = 0; i < n; i++) {
            findCluster(vecArray[i], clusterArray,k);
        }
        changed = reCalcMeans(clusterArray,k); /*re-calculate means. returns 1 if mean changed, 0 otherwise*/
        refreshClusters(clusterArray,k,k); /*init sum and count to zero's:*/
        iterCount++;
    }
}
/*functions*/
/*input: Array of clusters and it's length.number of coordinates.
output: None. All clusters sum and count is initalized to zero*/
void refreshClusters(Cluster** clusterArray,int k,int d) {
    int i, j;
    for (i = 0; i < k; i++) {
        clusterArray[i]->count = 0;
        /**loop through array of sum and init all to zero (0)**/
        for (j = 0; j < d; j++) {
            (clusterArray[i]->sum)[j] = 0;
        }
    }
}
/*input: 2 arrays(vectors coordinates)
output: squared distance between them*/
double distance(double *x, double *y, int d) {
    double sum = 0;
    int i;
    for (i = 0; i < d; i++) {
        sum += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sum;
}


/*deallocate memory*/
void freeMem(double** vecArray, Cluster** clusterArray,int n,int k) {
    int i;
    for(i=0;i<k;i++){
        free(clusterArray[i]->mean);
        free(clusterArray[i]->sum);
    }
    free(clusterArray);
    for(i=0;i<n;i++){
        free(vecArray[i]);
    }
    free(vecArray);
}

