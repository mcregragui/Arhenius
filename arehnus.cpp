#include <iostream>
#include<fstream>
#include<iomanip>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <string>
#include "arehnus.h"

using namespace std;
using namespace Eigen;



//fonction
double u(double x,double y, double t){
    return x*(1-x)*y*(1-y);
    //return sin(x)+cos(y);
}
double g(double x,double y, double t){
    return 0;
    //return cos(x)*x*x+cos(y);
}
double h(double x,double y, double t){
    if (t<=50){
        return -10000*t;
    }else{
        return -500000+9000*(t-50);
    }
    //return sin(x)-sin(y);
}
double g2(double x,double y, double t){
    //return 2*(1-2*x);
    return 0;
}
double h2(double x,double y, double t){
    //return 2*(1-2*y);
    return 0;
}
double f(double x,double y, double t){
    return 0;
    //return ;
    //return sin(x)+cos(y);
}


//systeme
int coefficient(int i,int j,int Nx){
    return (j-1)*Nx+i-1;
}
void remplissage(SparseMatrix<double>& matrice,VectorXd& b,double t,double deltax,double deltay,int Nx,int Ny){
    int coef;
    double x;
    double y;
    for (int j=1; j<=Ny;j++){
        for(int i=1; i<=Nx;i++){
            coef=coefficient(i,j,Nx);
            x=i*deltax;
            y=j*deltay;
            if(i==1){
                matrice.coeffRef(coef,coef+1)-=1./(deltax*deltax*1.);
                matrice.coeffRef(coef,coef)+=1./(deltax*deltax*1.);
                b(coef)+=f(x,y,t)-g(0,y,t)/(deltax*1.)-g2(0,y,t)/2.;
            }
            if(i==Nx){
                matrice.coeffRef(coef,coef)+=1./(deltax*deltax*1.);
                matrice.coeffRef(coef,coef-1)-=1./(deltax*deltax*1.);
                b(coef)+=f(x,y,t)+g(1,y,t)/(deltax*1.)-g2(1,y,t)/2.;//1=Lx
            }
            if(i!=Nx&&i!=1){
                matrice.coeffRef(coef,coef)+=2./(deltax*deltax*1.);
                matrice.coeffRef(coef,coef-1)-=1./(deltax*deltax*1.);
                matrice.coeffRef(coef,coef+1)-=1./(deltax*deltax*1.);
                b(coef)+=f(x,y,t);
            }


            if(j==1){
                matrice.coeffRef(coef,coef)+=1./(deltay*deltay*1.);
                matrice.coeffRef(coef,coef+Nx)-=1./(deltay*deltay*1.);
                b(coef)-=h(x,0,t)/(deltay*1.)+h2(x,0,t)/2.;
            }
            if(j==Ny){
                matrice.coeffRef(coef,coef)+=1./(deltay*deltay*1.);
                matrice.coeffRef(coef,coef-Nx)-=1./(deltay*deltay*1.);
            }
            if(j!=Ny&&j!=1){
                matrice.coeffRef(coef,coef)+=2./(deltay*deltay*1.);
                matrice.coeffRef(coef,coef-Nx)-=1./(deltay*deltay*1.);
                matrice.coeffRef(coef,coef+Nx)-=1./(deltay*deltay*1.);
            }
        }
    }
}
void remplissagep(SparseMatrix<double>& matrice,VectorXd& b,VectorXd rhon, VectorXd rhonp,double t,double deltax,double deltay,int Nx,int Ny){
   
    int coef;
    double x;
    double y;
    VectorXd s(Nx*Ny);
    s=vmpointg(rhon,rhonp,Nx,Ny, t);
    for (int j=1; j<=Ny;j++){
        for(int i=1; i<=Nx;i++){
            coef=coefficient(i,j,Nx);
            x=i*deltax;
            y=j*deltay;
          
            if(i==1){
                matrice.coeffRef(coef,coef+1)-=1./(deltax*deltax*1.);
                matrice.coeffRef(coef,coef)+=1./(deltax*deltax*1.)+s(coef);
                b(coef)+=f(x,y,t)-g(0,y,t)/(deltax*1.)-g2(0,y,t)/2.;
            }
            if(i==Nx){
                matrice.coeffRef(coef,coef)+=1./(deltax*deltax*1.)+s(coef);
                matrice.coeffRef(coef,coef-1)-=1./(deltax*deltax*1.);
                b(coef)+=f(x,y,t)+g(1,y,t)/(deltax*1.)-g2(1,y,t)/2.;//1=Lx
            }
            if(i!=Nx&&i!=1){
                matrice.coeffRef(coef,coef)+=2./(deltax*deltax*1.)+s(coef);
                matrice.coeffRef(coef,coef-1)-=1./(deltax*deltax*1.);
                matrice.coeffRef(coef,coef+1)-=1./(deltax*deltax*1.);
                b(coef)+=f(x,y,t);
            }


            if(j==1){
                matrice.coeffRef(coef,coef)+=1./(deltay*deltay*1.)+s(coef);
                matrice.coeffRef(coef,coef+Nx)-=1./(deltay*deltay*1.)-s(coef+Nx);
                b(coef)-=h(x,0,t)/(deltay*1.)+h2(x,0,t)/2.;
            }
            if(j==Ny){
                matrice.coeffRef(coef,coef)+=1./(deltay*deltay*1.)+s(coef);
                matrice.coeffRef(coef,coef-Nx)-=1./(deltay*deltay*1.);
            }
            if(j!=Ny&&j!=1){
                matrice.coeffRef(coef,coef)+=2./(deltay*deltay*1.)+s(coef);
                matrice.coeffRef(coef,coef-Nx)-=1./(deltay*deltay*1.);
                matrice.coeffRef(coef,coef+Nx)-=1./(deltay*deltay*1.)-s(coef+Nx);
            }
        }
    }
}
void initial_valeur(VectorXd& T,int Nx,int Ny,double valeur){
    for(int i=0;i<Nx*Ny;i++){
        T(i)=valeur;
    }
}
void remplissage_rho(SparseMatrix<double>& matrice_rho,VectorXd& b_rho,VectorXd T,double t,double deltax,double deltay,int Nx,int Ny){
    int coef;
    double x;
    double y;
    double rho_v=1500.;
    double rho_p=1000.;
    double A=1000.;
    double Ta=6000.;
    for (int j=1; j<=Ny;j++){
        for(int i=1; i<=Nx;i++){
            coef=coefficient(i,j,Nx);
            matrice_rho.coeffRef(coef,coef)=(rho_v*A*exp(-Ta/(T(coef)*1.)))/((rho_v-rho_p)*1.);
            b_rho(coef)=(rho_v*rho_p*A*exp(-Ta/(T(coef)*1.)))/((rho_v-rho_p)*1.);
        }
    }
}
void rho_fois_cp(VectorXd& rho_cp,VectorXd rho_etoile,int Nx,int Ny){
    double C_p=1500.;
    double C_v=1000.;
    double rho_p=1000.;
    double rho_v=1500.;
    double x;
    for(int i=0;i<Nx*Ny;i++){
        x=(rho_v-rho_etoile(i))/((rho_v-rho_p)*1.);
        rho_cp(i)=(1-x)*rho_v*C_v+x*C_p*rho_p;
    }
}
void v_m(Eigen::SparseMatrix<double>& matrice_rho_cp,Eigen::VectorXd rho_cp,int Nx,int Ny){
    for(int i=0;i<Nx*Ny;i++){
        matrice_rho_cp.coeffRef(i,i)=rho_cp(i);
    }
}
double mpointg(int i , int j ,Eigen::VectorXd rhon,Eigen::VectorXd rhonp,int Nx,int Ny,double t)
{
    /*int n=Nx*Ny;
    rhon.resize(n);
    rhonp.resize(n);*/
    double g(0.),m1(0.),m(0.);
    double coef;
    for (int k=1; k<=j;k++)
    {
        
        coef=coefficient(i,k,Nx);
        g+=(rhonp(coef)-rhon(coef)); 
        
    }
    //cout<< "here "<<g<<endl;
    double coef0=coefficient(i,1,Nx);
    double coef1=coefficient(i,j,Nx);
    //cout<<"here "<<((rhonp(coef0)-rhon(coef0))+(rhonp(coef1)-rhon(coef1)))/2.<<" g: "<<g<<" res: "<<((rhonp(coef0)-rhon(coef0))+(rhonp(coef1)-rhon(coef1)))/2.-g<<endl;
    m1=-g;
    cout<<m1<<endl;
   // m=m1;
    
    return m1;
}
Eigen::VectorXd vmpointg(Eigen::VectorXd rhon,Eigen::VectorXd rhonp,int Nx,int Ny,double t)
{
    double coef;
    Eigen::VectorXd mp(Nx*Ny);
    for (int j = 1; j <=Ny; j++)
    {
        for (int i = 1; i <= Nx; i++)
        {
            coef=coefficient(i,j,Nx);
            cout<<"here4"<<endl;
            mp(coef)=1000000.*mpointg(i,j,rhon,rhonp, Nx, Ny,t);
            cout<<"here5"<<endl;
        }
        return mp;
    }
    
}
//outil
void print_x(VectorXd v,string nom,int j,int Nx,int Ny,double deltax,double t){
    ofstream mon_flux;
    string name_file=nom+" la trache "+to_string(j) +"  l'instant "+to_string(t);
    mon_flux.open(name_file, ios::out);
    double x;
    int coef;
    for(int i=1;i<=Nx;i++){
        x=i*deltax;
        coef=coefficient(i,j,Nx);
        mon_flux<<x<<" "<<std::setprecision(16)<< v(coef) << endl;
    }
    mon_flux.close();
}
void print_y(VectorXd v,string nom,int i,int Nx,int Ny,double deltay,double t){
    ofstream mon_flux;
    string name_file=nom+" la trache "+to_string(i) +"l'instant "+to_string(t);
    mon_flux.open(name_file, ios::out);
    double y;
    int coef;
    for(int j=1;j<=Ny;j++){
        y=j*deltay;
        coef=coefficient(i,j,Nx);
        mon_flux<<y<<" "<<std::setprecision(16)<< v(coef) << endl;
    }
    mon_flux.close();
}
