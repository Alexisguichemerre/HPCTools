#include <iostream>
#include <time.h> 
#include <math.h>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace std;


// Function :  Decomposition of the matrix a with the LU method without pivoting.
// Input : A matrixXd a , an integer n
// Output : A vector of matrixXd

vector<MatrixXd> lu(MatrixXd a, int n)
{
    MatrixXd l = MatrixXd::Identity(n,n);
    MatrixXd u = MatrixXd::Zero(n,n);
	vector<MatrixXd > vec(2);
    int i = 0, j = 0, k = 0;
    int p ;
    for (i = 0; i < n; i++)
    {
        p =i+1;
        for (j = 0; j < p ; j++)
        {
            u(j,i) = a(j,i);
            for (k = 0; k < j; k++)
            {
                u(j,i) = u(j,i) - (l(j,k) * u(k,i)) ;
            }
        }
        for (j = i; j < n ; j++)
        {

            l(j,i) = a(j,i);
            for (k = 0; k < i; k++)
            {
                l(j,i) = l(j,i) - (l(j,k) * u(k,i));
            }
            l(j,i)=l(j,i)/u(i,i);
            
        }
    }

    vec[0]=l;
    vec[1]=u;
    return vec;
}




// Function : Resolution of the system LZ=B
//Input : A matrixXd, A matrixXd, an integer
//Output : A matrixXd

MatrixXd Descent(MatrixXd L, MatrixXd B, int n){
    MatrixXd X(n,n);
    double val;
   for (int i = 0; i < n; i++)
   {
       for (int j = 0; j < n; j++)
          {
            val = 0.0;
            for (int k = 0; k < j; k++)
            {
                val+=L(j,k)*X(k,i);
            }
              X(j,i)=(B(j,i)-val)/L(j,j);
          }  
    }
return X;
}




// Function : Resolution of the system UX=Z
//Input : A matrixXd, A matrixXd, an integer
//Output : A matrixXd
MatrixXd Mount(MatrixXd U, MatrixXd Z, int n){
    MatrixXd X(n,n);
    int p = n-1;
    double val;
   for (int i = 0; i < n; i++)
   {
       for (int j = p; j >=0; j--)
          {
            val = 0.0;
            for (int k = j+1; k < n; k++)
            {
                val+=U(j,k)*X(k,i);
            }
              X(j,i)=(Z(j,i)-val)/U(j,j);
          }  
    }
return X;
}


// Function : Resolution of the initial system of equation AX=B
// Input : A MatrixXd, A MatrixXd
// Output : A matrixXd

MatrixXd Resolve(MatrixXd A,MatrixXd B){
    int Nb=A.cols();
    cout << Nb << endl;
    vector<MatrixXd> vec =lu(A,Nb);
    MatrixXd L;
    MatrixXd U;
    L=vec[0];
    U=vec[1];

    MatrixXd S = L*U;
    MatrixXd res;

    res = Descent(L,B,Nb);
    MatrixXd resf;

    resf = Mount(U,res,Nb);

return resf;    

}



int main(int argi, char** argc){

MatrixXd A=MatrixXd::Random(10,10);
MatrixXd B=MatrixXd::Random(10,10);

MatrixXd test = Resolve(A,B);

cout << test << endl;


}