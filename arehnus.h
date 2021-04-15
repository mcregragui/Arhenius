#include "Eigen/Dense"
#include "Eigen/Sparse"



double u(double x,double y, double t);
double g(double x,double y, double t);
double h(double x,double y, double t);
double f(double x,double y, double t);


int coefficient(int i,int j,int Nx);
void remplissage(Eigen::SparseMatrix<double>& matrice,Eigen::VectorXd& b,double t,double deltax,double deltay,int Nx,int Ny);
void remplissage_rho(Eigen::SparseMatrix<double>& matrice_rho,Eigen::VectorXd& b_rho,Eigen::VectorXd T,double t,double deltax,double deltay,int Nx,int Ny);
void initial_valeur(Eigen::VectorXd& T,int Nx,int Ny,double valeur);
void rho_fois_cp(Eigen::VectorXd& rho_cp,Eigen::VectorXd rho_etoile,int Nx,int Ny);
void v_m(Eigen::SparseMatrix<double>& matrice_rho_cp,Eigen::VectorXd rho_cp,int Nx,int Ny);
double mpointg(int , int ,Eigen::SparseMatrix<double>,int Nx,int Ny,double t);
Eigen::VectorXd vmpointg(Eigen::VectorXd rhon,Eigen::VectorXd rhonp,int Nx,int Ny,double t);
void remplissagep(Eigen::SparseMatrix<double>& matrice,Eigen::VectorXd& b,Eigen::VectorXd rhon, Eigen::VectorXd rhonp,double t,double deltax,double deltay,int Nx,int Ny);
void print_x(Eigen::VectorXd v,std::string nom,int j,int Nx,int Ny,double deltax,double t);
void print_y(Eigen::VectorXd v,std::string nom,int i,int Nx,int Ny,double deltay,double t);