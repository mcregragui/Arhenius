#include <iostream>
#include<fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseQR"
#include "Eigen/SparseCholesky"
#include <string>
#include "arehnus.h"

using namespace std;
using namespace Eigen;





int main()
{
    cout<<"bienvenu dans le code de calcul pyrolyse"<<endl;
    int Nx=5;
    double Lx=0.0025;
    double deltax=Lx/((Nx+1)*1.);
    
    
    int Ny=250;
    double Ly=0.01;
    double deltay=Ly/((Ny+1)*1.);

    double tf=100;
    double dt=0.1;
    double t=dt;
    double x,y;
    int coef;
    int nb=int(ceil(tf/dt));

    //valeur
    double T0=293.;
    double rho0=1500.;
    double Lm=3000000;

   
    VectorXd T(Nx*Ny);
    VectorXd rho(Nx*Ny);
    VectorXd rho1(Nx*Ny);
    VectorXd rho_etoile(Nx*Ny);
    VectorXd rho_cp(Nx*Ny);
    VectorXd rho_cp_1(Nx*Ny);
    VectorXd T1(Nx*Ny);

    SparseMatrix<double> matrice(Nx*Ny,Nx*Ny),id(Nx*Ny,Nx*Ny);
    SparseMatrix<double> matrice_rho(Nx*Ny,Nx*Ny);
    SparseMatrix<double> matrice_rho_cp(Nx*Ny,Nx*Ny);
    SparseMatrix<double> matrice_rho_cp_1(Nx*Ny,Nx*Ny);
    VectorXd b(Nx*Ny);
    VectorXd b_rho(Nx*Ny);
    matrice.setZero();
    b.setZero();
    


    
    initial_valeur(T,Nx,Ny,T0);
    initial_valeur(rho,Nx,Ny,rho0);
    id.setIdentity();
    matrice_rho_cp.setZero();
    VectorXd temp(Nx*Ny);
    VectorXd temp1(Nx*Ny);
    temp1=rho;
    for(int n=1;n<=nb;n++){

        
        //prediction
        
            remplissage_rho(matrice_rho,b_rho,T,t,deltax,deltay,Nx,Ny);
            SimplicialLLT <SparseMatrix<double> > solver_rho;
            solver_rho.compute(dt*matrice_rho+id);
            rho_etoile=solver_rho.solve(rho+dt*b_rho);
            rho_fois_cp(rho_cp,rho,Nx,Ny);
            rho_fois_cp(rho_cp_1,rho_etoile,Nx,Ny);
            v_m(matrice_rho_cp_1,rho_cp_1,Nx,Ny);
            v_m(matrice_rho_cp,rho_cp,Nx,Ny);

        
        //calcul T
        b.setZero();
        matrice.setZero();
        if(n==1)
        {
            remplissage(matrice,b,t,deltax,deltay,Nx,Ny);
        }
        else
        {
            remplissagep(matrice,b,temp,temp1,t,deltax,deltay,Nx,Ny);
        }
        b+=dt*Lm*(matrice_rho*rho_etoile-b_rho);
        SimplicialLLT <SparseMatrix<double> > solver;
        solver.compute(dt*matrice+matrice_rho_cp_1);
        T1=solver.solve(matrice_rho_cp*T+dt*b+T0*(rho_cp_1-rho_cp));

        if(n==1){
            //cout<<matrice_rho_cp<<endl;
            cout<<T1<<endl;
        }



        //correction
        remplissage_rho(matrice_rho,b_rho,T1,t,deltax,deltay,Nx,Ny);
        SimplicialLLT <SparseMatrix<double> > solver_rho1;
        solver_rho1.compute(dt*matrice_rho+id);
        rho1=solver_rho1.solve(rho+dt*b_rho);
        if(n==1)
        {
          temp1=rho1;  
        }
        else{
        temp=temp1;
        temp1=rho1;
        }
        if(n%100==0){
            print_y(T1,"temps",2,Nx,Ny,deltay,t);
            print_y(rho1,"rho",2,Nx,Ny,deltay,t);
        }
        T=T1;
        rho=rho1;
        t+=dt;
    }
   
    
    return 0;
}