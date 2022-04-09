#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <ctime>
#include <math.h>


using namespace std;


///////////////////////////////////// Initialisation /////////////////////////////////

double *zero(int N)       //vecteur 0
{
    double *b = new double [N];
    for (int i=0;i<N;i++){
        b[i]=0;
    }
    return b;
}

double *one(int N)     //vecteur qui ne contient que des 1 
{
    double *b = new double [N];
    for (int i=0;i<N;i++){
        b[i]=1;
    }
    return b;
}


double *linspace(double a, double b, int N)  // discretisation de l'intervalle
{
	double dx = (b - a) / (N - 1);
	double *v = zero(N);
	for (int i = 0; i < N; i++) {
		v[i] = a+i*dx;
	}
	return v;
}

double **Zero(int N, int M) //créer une matrice de taille N*M remplie avec des zéros
{
    double **A = new double *[N];
    for (int i=0;i<N;i++){
        A[i] = new double [M];
        for(int j=0;j<M;j++){
            A[i][j]=0;
            
        }
    }
    return A;
}

double ps(double*u,double*v,int N)     //produit scalaire
{
    double r=0;
    for(int i=0;i<N;i++){
        r+=u[i]*v[i];
    }
    return r;
}


double **pm(double **A,double **B,int N)    //produit matriciel
{
    double **C = Zero(N,N);
    for (int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for (int k=0;k<N;k++){
                C[i][j] += A[i][k]*B[k][j];
            }
            
        }
    }
    return C;
}

double *pmv(double **A,double *v,int N)   // produit matrice vecteur
{
    double *u = zero(N);
    for (int i=0;i<N;i++){
        u[i] = ps(A[i],v,N);
    }
    return u;
}

double **am(double **A, double **B,int N)  // addition matricelle
{
    double **C = Zero(N,N);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            C[i][j] = A[i][j]+B[i][j];
        }
    }
    return C;
}

double *av(double *u,double *v,int N)    // addition vectorielle
{
    double *r = zero(N);
    for (int i=0;i<N;i++){
        r[i] = u[i]+v[i]; 
    }
    return r;
}

double *dv(double *u,double *v,int N)    // difference entre les 2 vecteurs
{
    double *r = zero(N);
    for (int i=0;i<N;i++){
        r[i] = u[i]-v[i]; 
    }
    return r;
}

double **prm(double **A, double b, int N)   //produit reel matrice
{
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            A[i][j] *= b;
        }
    }
    return A;
}

double *prv(double *v,double a,int N)    //produit reel vecteur
{
    for (int i=0;i<N;i++){
        v[i] *= a;
    }
    return v;
}

double norme(double *u, int N)   // norme 2
{
    double r = 0;
    for(int i=0;i<N;i++){
        r += u[i]*u[i];
    }
    r = sqrt(r);
    return r;
}


double max(double *u,int N) //element max du vecteur
{
    double max = -999;
    for (int i=0;i<N;i++){
        if (u[i]>max){
            max = u[i];
        }
    }
    return max;
}


void print_m(double**A,int N,int M)    //affichage de matrice
{
    for( int i =0; i<N ; i++){
        cout<<"ligne"<<i<<": " ;
        for ( int j =0; j<M ; j++){
            cout<<A[i][j]<<" " ;
        }
    cout<<endl ;
    }
}


void print_v(double*A,int N)    //affichage de vecteur
{
    cout<<"[";
    for( int i =0; i<N-1 ; i++){
        cout<<A[i]<<"," ;
    }
    cout<<A[N-1]<<"]";
}



double *jacobi(double **A,double *b, double *x0,int N)    //methode de jacobi
{
    int maxiter = 10000;
    double tol = 0.00000001;
    double **S = Zero(N,N);
    double **T = Zero(N,N);
    double *x = x0;
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            if (i==j){
                S[i][j] = A[i][j];
            }
            else{
                T[i][j] = A[i][j]*(-1);
            }
        }
    }
    double **invS = Zero(N,N);
    for(int i=0;i<N;i++){
        invS[i][i] = 1/S[i][i];
    }
    double *xx = zero(N);
    for (int i=0;i<maxiter;i++){
        xx = av(pmv(pm(invS,T,N),x,N),pmv(invS,b,N),N);
        if ((norme(av(xx,prv(x,-1,N),N),N))-tol < 0){
            break;
        }
        x = prv(xx,1,N);
        if (i == maxiter-1){
            return zero(N);
        }
    }
    double *sol = xx;
    return sol;
}



/////////////////////////////////////// EXO 1 /////////////////////////////////////////

double f(double x)     //fonction de teste
{
    double r = 4*pow(x,3)+ 3*x*x + 2*x+1;
    return r;
}

double *moment(double a,double b,int N,double (*f)(double a)){
    double *x = linspace(a,b,N);
    double *M = one(N-2);    //valeur initiale
    double h = (b-a)/(N-1);  //pas
    double *y = zero(N-2);   //vecteur pour stocker la solution exacte
    double *yy = zero(N);  //solution exacte avec le bord
    yy[0] = (f(a));
    yy[N-1] = (f(b));
    for (int i=1;i<N-1;i++){
        y[i-1] = f(x[i]);
    }
    double  **A = Zero(N-2,N-2);   //matrice tridiagonale
    for (int i=0;i<N-2;i++){
        for (int j=0;j<N-2;j++){
            if (i==j){
                A[i][j] = 4;
            }
            if(i-j==1||i-j==-1){
                A[i][j] = 1;
            }
        }
    }
    double **B = Zero(N-2,N-2);  //matrice tridiagonale
    for (int i=0;i<N-2;i++){
        for (int j=0;j<N-2;j++){
            if (i==j){
                B[i][j] = -12/h/h;
            }
            if(i-j==1||i-j==-1){
                B[i][j] = 6/h/h;
            }
        }
    }
    double *u = pmv(B,y,N-2);  
    double *R = jacobi(A,u,M,N-2); //methode d'inversion pour Ax = u
    double *sol = zero(N);   //vecteur pour stocker la solution
    for (int i=1;i<N-1;i++){
        sol[i] = R[i-1];
    }
    return sol;
}


double interpol(double *M,double a,double b,int N,double x ,double (*f)(double a)){
    double *d = linspace(a,b,N);
    int r = 0;
    for (int i=0;i<N;i++){
        if (x<d[i]){
            r = 1*i;    //trouver la position de x dans l'intervalle
            break;
        }
    }
    double h = (b-a)/(N-1);  //pas
    double Sf = (-M[r-1]/6/h)*pow(x-d[r],3)+M[r]/h*pow(x-d[r-1],3)+
    ((f(d[r])-f(d[r-1]))/h+h/6*(M[r-1]-M[r]))*(x-d[r-1])+f(d[r-1])-M[r-1]*h*h/6;  //Sf
    return Sf;
}

//////////////////////////////////// EXO 2 /////////////////////////////////////

double *euler_explicite(double T, double N){
    double *y = zero(N+1);
    y[0]=1;
    double h = T/N;
    for (int i=0;i<N;i++){
        double ti=i*h;
        y[i+1]=(1-2*h*ti)*y[i];
    }
    return y;
}
double *euler_implicite(double T, double N){
    double *y =zero(N+1);
    double h = T/N;
    y[0]=1;
    for (int i=0;i<N;i++){
        double ti1=(i+1)*h;
        y[i+1]=y[i]* (1/(1+2*h*ti1));
    }
    return y;
}
double *crank_nicholson(double T,double N){
    double h=T/N;
    double *y =zero(N+1);
    y[0]=1;
    for (int i=0;i<N;i++){
        double ti=i*h;
        double ti1=(i+1)*h;
        y[i+1]=y[i]*(1+h*ti)/(1+h*ti1);
    }
    return y;

}
double *Heun(double T,double N){
    double h=T/N;
    double *y =zero(N+1);
    y[0]=1;
    for (int i=0;i<N;i++){
        double ti=i*h;
        double ti1=(i+1)*h;
        y[i+1]=y[i]*(1-(ti+ti1)*h)+2*h*h*ti*ti1;
    }
    return y;

}

/////////////////////////////// EXO 3 /////////////////////////////////

double*F(double*V,double c){
    // c = masse de la terre * Constante de gravitation universelle
    double*U=zero(4);
    U[0]=V[2];
    U[1]=V[3];
    U[2]=c*V[0]/pow((V[0]*V[0]+V[1]*V[1]),1.5);
    U[3]=c*V[1]/pow((V[0]*V[0]+V[1]*V[1]),1.5);
    return U;
}


double trajectoire(double d, double v_init, double N,double*(*F)(double*V,double c),double c){
    // d est la distance initiale 6500km
    // v_init est la vitesse initiale du satellite (en km)
    // 
    double h=8600/N; // temps
    double*U=zero(4);
    U[0]=d;
    U[3]=v_init;
    double* dis = zero(N);  //vecteur pour stocker la distance
    
    for (int i=0;i<N;i++){
        double*FUi=F(U,c); //F(U_i)
        double*hFUi=prv(FUi,h,4); //h*F(U_i)
        double*FUih=F((av(U,hFUi,4)),c); //F(U_i+h*F(U_i))
        U=av(U,prv(av(FUi,FUih,4),h/2,4),4);
        dis[i] = sqrt(U[0]*U[0]+U[1]*U[1]);
    }
    return max(dis,N);
}

/////////////////////////////// EXO 4 Partie I/////////////////////////////////

double c(double t){ //fonction c(t)=t^2+2
    return t*t+2;
}

double f2(double t){ //fonction f(t)=t^6+5t^4-5t^2-4,
    return pow(t,6)+5*pow(t,4)-5*t*t-4;
}
double* resouScd(double N, double alpha, double beta, double (*c)(double x), 
double(*f2)(double y)){
    double h =1/N;//discretisation de l'intervalle [0,1]
    double* U=zero(N+1);
    U[0]=alpha;
    U[1]=alpha+beta*h; //Valeur initiale
    for (int i=1;i<N;i++){
        double t_i=i*h;
        U[i+1]=h*h*(c(t_i)*U[i]-f2(t_i))+2*U[i]-U[i-1];
    }
    return U;
}

//fonction u(t)=t^4+3*t^2+1, u"(t)=12t^2+6, alpha=u(0)=1, beta=u'(0)=0
double* fonction_U(double N){
    double h=1/N;
    double*U=zero(N+1);
    U[0]=1;
    for(int i=1;i<N+1;i++){
        double t_i=i*h;
        U[i]=pow(t_i,4)+3*t_i*t_i+1;
    }
    return U;
}

/////////////////////////////// EXO 4 Partie II /////////////////////////////////

double f3(double t){ //fonction f(t)=t^4-t^3+2t^2-2t-2,
    return pow(t,4)-pow(t,3)+2*t*t-2*t-2;
}

double* resouScd2(double N, double (*c)(double x), double(*f3)(double y)){
    double h =1/N;//discretisation de l'intervalle [0,1]
    double* U=zero(N+1);   //stocker la solution du bord
    double **A = Zero(N-1,N-1);  
    double *b = zero(N-1);
    for (int i=0;i<N-1;i++){
        if (i==0){
            double t_i=(i+1)*h;
            A[i][i] = 2/h/h+c(t_i);
            A[i][i+1] = -1/h/h;
            b[i] = f3(t_i);
        }
        if (i>0 && i<N-2){
            double t_i=(i+1)*h;
            A[i][i] = 2/h/h+c(t_i);
            A[i][i+1] = -1/h/h;
            A[i][i-1] = -1/h/h;
            b[i] = f3(t_i);
        }
        if (i==N-2){
            double t_i=(i+1)*h;
            A[i][i] = 2/h/h+c(t_i);
            A[i][i-1] = -1/h/h;
            b[i] = f3(t_i);
        }
    }
    double *x0 = one(N-1);
    double *uu = jacobi(A,b,x0,N-1);  // la methode de jacobi
    for(int i=0;i<N-1;i++){
        U[i+1] = uu[i];
    }
    return U;
}
 
//fonction u(t)=t^2-t, u(0)=0, u(1)=0
double* fonction_U2(double N){
    double h=1/N;
    double*U=zero(N+1);
    U[0]=0;
    for(int i=1;i<N+1;i++){
        double t_i=i*h;
        U[i]=pow(t_i,2)-t_i;
    }
    return U;
}






int main(){


/////////////// EXO 1 //////////////
   /*
   int a = -5;
   int b = 5;
   int N = 100;
   double *sol = moment(a,b,N,&f);
   double *s = new double[4];
   s[0] = -4.3; s[1] = 0; s[2] = 3; s[3] = 4.99;
   for (int i=0;i<4;i++){
       double Sf = interpol(sol,a,b,N,s[i],&f);
        cout<<"Solution approchee:"<<Sf<<endl;
        cout<<"Solution exacte:"<<f(s[i])<<endl;;
   }
   */

  //////////// EXO 2 ///////////////
    /*
    double T = 3 ;
    double N = 50 ;
    print_v ( euler_explicite (T , N ) , N +1) ;
    cout<<""<<endl;
    print_v ( euler_implicite (T , N ) , N +1) ;
    cout<<""<<endl;
    print_v ( Heun (T , N ) , N +1) ;
    cout<<""<<endl;
    print_v ( crank_nicholson (T , N ) , N +1);
    */

  //////////// EXO 3 //////////////

  /*
    double d=6500;
    double v1=8;
    double v2 = 10;
    double v3 = 12;
    double N=10000;
    double c=-5.9722*6.6743*pow(10,4);
    double traj1 = trajectoire(d,v1,N,&F,c);
    double traj2 = trajectoire(d,v2,N,&F,c);
    double traj3 = trajectoire(d,v3,N,&F,c);
    cout<<"apogee de la trajectoire est : "<<traj1<<endl;
    cout<<"apogee de la trajectoire est : "<<traj2<<endl;
    cout<<"apogee de la trajectoire est : "<<traj3<<endl;
  */


  //////////// EXO 4  Partie I//////////////

  /*
    double N = 500 ;
    double alpha = 1 ;
    double beta = 0 ;
    print_v ( resouScd (N , alpha , beta ,& c ,& f2 ) , N +1) ;
    cout<<""<<endl;
    print_v (fonction_U(N), N+1);
    cout<<""<<endl;
    print_v ( dv(fonction_U ( N ),resouScd (N , alpha , beta ,& c ,& f2 ),N+1) , N+1);
  */


  //////////// EXO 4  Partie II //////////////
  
  /*
    double N = 50;
    double *u = resouScd2(N,&c,&f3);  //solution approchee
    print_v(u,N+1);
    cout<<""<<endl;
    print_v(fonction_U2(N), N+1);  //solution exacte
    cout<<""<<endl;
    print_v(dv(u,fonction_U2(N),N+1) , N+1); //visualisation d'erreur
  */

    return 0;
}