#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h> 
#include <math.h>
#include <vector>
// #include "mkl_lapacke.h"

using namespace std;

double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}





void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", name);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}



int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if (bref[i]!=b[i]) return 0;
    }
    return 1;
}

//Function : Initialisation of a " matrix  with pointer"
// Input : A int size
// Ouput : A double** correspond to the matrix
double** init(int size) {
    double** a = new double*[size];
    for (int i = 0; i < size; ++i)

        {
            a[i]= new double[size];
        }
    return a;
}


//Function : Initialisation of a " zero matrix"
// Input : A int size
// Ouput : A double** correspond to the zero matrix
double** zero(int size) {
    double** a = new double*[size];
    for (int i = 0; i < size; ++i)

        {
            a[i]= new double[size];
        }
    for (int i = 0; i < size; ++i)
        {
        for (int j = 0; j <size; ++j)
            {
                a[i][j]= 0.0; 
            }
        }
    return a;
}


double** matrixgenerate(double *matrix, int size) {
    double** C=init(size);
    int p =0;
    for (int i = 0; i < size ; i++)
    {
        p=i*size;
        for (int j = 0; j < size ; j++)
        {
            C[i][j]=matrix[p+j];
        }
    }
    return C;
}

//Function : Initialisation of a " identity matrix"
// Input : A int size
// Ouput : A double** correspond to the identity matrix
double** Ident(int size) {
    double** a = new double*[size];
    for (int i = 0; i < size; ++i)

        {
            a[i]= new double[size];
        }
    for (int i = 0; i < size; ++i)
        {
        for (int j = 0; j <size; ++j)
            {
                a[i][j]= 0.0; 
            }
        }
    return a;
}


// Permit to  display the matrix
void disp(double** A, int size){
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << A[i][j] << " " ;
        }
        cout << endl;
    }
}



// Function :  Decomposition of the matrix a with the LU method without pivoting.
// Input : A double** a , an integer n
// Output : A vector of double**

vector<double** > lu(double** a, int n)
{
    double** l = Ident(n);
    double** u = zero(n);
    vector<double** > vec(2);
    int i = 0, j = 0, k = 0;
    int p ;
    for (i = 0; i < n; ++i)
    {
        p =i+1;
        for (j = 0; j < p ; ++j)
        {
            u[j][i] = a[j][i];
            for (k = 0; k < j; ++k)
            {
                u[j][i] = u[j][i] - (l[j][k] * u[k][i]) ;
            }
        }
        for (j = i; j < n ; ++j)
        {

            l[j][i] = a[j][i];
            for (k = 0; k < i; ++k)
            {
                l[j][i] = l[j][i] - (l[j][k] * u[k][i]);
            }
            l[j][i]=l[j][i]/u[i][i];
            
        }
    }

    vec[0]=l;
    vec[1]=u;
    return vec;
}


// Function : Resolution of the system LZ=B
//Input : A double**, A double**, an integer
//Output : A double**

double** Descent(double** L, double** B, int n){
    double** X = init(n);
    double val;
   for (int i = 0; i < n; ++i)
   {
       for (int j = 0; j < n; ++j)
          {
            val = 0.0;
            for (int k = 0; k < j; ++k)
            {
                val+=L[j][k]*X[k][i];
            }
              X[j][i]=(B[j][i]-val)/L[j][j];
          }  
    }
return X;
}


// Function : Resolution of the system UX=Z
//Input : A double**, A double**, an integer
//Output : A double**
double** Mount(double** U, double** Z, int n){
    double** X=init(n);
    int p = n-1;
    double val;
   for (int i = 0; i < n; ++i)
   {
       for (int j = p; j >=0; --j)
          {
            val = 0.0;
            for (int k = j+1; k < n; ++k)
            {
                val+=U[j][k]*X[k][i];
            }
              X[j][i]=(Z[j][i]-val)/U[j][j];
          }  
    }
return X;
}



// Function : Resolution of the initial system of equation AX=B
// Input : A MatrixXd, A MatrixXd
// Output : A matrixXd

double** Resolve(double** A,double** B,int n){
    vector<double**> vec =lu(A,n);
    double** L;
    double** U;
    L=vec[0];
    U=vec[1];

    double**  res;

    res = Descent(L,B,n);
    double** resf;

    resf = Mount(U,res,n);

return resf;    

}

// Function : Product between 2 matrices.
// Input : A double** A, A double** C , int size
// Output : A double** B

double** prod(double** A, double** C,int size){
    double** B = zero(size);
    double compt=0.0;
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            for (int k = 0; k < size; ++k)
            {
                compt = compt + (A[i][k]*C[k][j]);
            }
            B[i][j]=compt;
            compt=0.0;
        }
    }
    return B;
}


double** randMat(int n){
    double** mat = init(n);
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            double val = (double)rand() / RAND_MAX;
            mat[i][j] = val;
        }
    }
    return mat;
}


    int main(int argc, char *argv[])
    {

        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size);
        aref = generate_matrix(size);        
        b = generate_matrix(size);
        bref = generate_matrix(size);


        // Generation of a 2 random matrix with the library eigen.
        // double** A=matrixgenerate(a,size);
        //  double** B=matrixgenerate(b,size);
        double** A=randMat(size);
        double** B=randMat(size);



        // print_matrix("A", a, size);
        // print_matrix("B", b, size);

        // Using MKL to solve the system
        // MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        // MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        // info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        // printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        tStart = clock();    
        // MKL_INT *ipiv2 = (MKL_INT *)malloc(sizeof(MKL_INT)*size);       
        double** test = Resolve(A,B,size);
        // disp(prod(A,test,size),size);
        // cout << "***" << endl;
        // disp(B,size);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        if (check_result(bref,b,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        
        // print_matrix("X", b, size);
        // print_matrix("Xref", bref, size);
        return(0);
    }
