#include "charge.h"
#include <iomanip>
#include "config.h"
#include "mmm/constants.h"
#include "mesh/Field.h"
#include "mesh/VectorField.h"
#include "mesh/RectangularMesh.h"
#include "Magneto.h"
#include "Benchmark.h"
#include <iostream>

double topology_continuous_helper(int i, int j, int k, bool px, bool py, int dim_x, int dim_y,VectorField::const_accessor &M_acc,Field::ro_accessor &Ms_acc){
    double Msij = 0;
    Vector3d Mij,Mx,My;
    Msij = Ms_acc.at(i,j,k);
    if(Msij == 0.0){
        return 0;
    }
    Mij = M_acc.get(i,j,k);

    if(j==0){
        if(py){
            My = (M_acc.get(i,j+1,k) - M_acc.get(i,dim_y-1,k))/2;
        }else{
            My = M_acc.get(i,1,k) - Mij;
        }
    }else if(j == dim_y-1){
        if(py){
            My = (M_acc.get(i,0,k) - M_acc.get(i,j-1,k))/2;
        }else{
            My = Mij - M_acc.get(i,dim_y-2,k);
        }
    }else{
        My = (M_acc.get(i,j+1,k) - M_acc.get(i,j-1,k))/2;
    }

    if(i==0){
        if(px){
            Mx = (M_acc.get(i+1,j,k) - M_acc.get(dim_x-1,j,k))/2;
        }else{
            Mx = M_acc.get(1,j,k) - Mij;    
        }
    }else if(i == dim_x-1){
        if(px){
            Mx = (M_acc.get(0,j,k) - M_acc.get(i-1,j,k))/2;
        }else{
            Mx = Mij - M_acc.get(dim_x-2,j,k);    
        }
    }else{
        Mx = (M_acc.get(i+1,j,k) - M_acc.get(i-1,j,k))/2;
    }
    return dot(Mij,cross(Mx,My))/(Msij*Msij*Msij);
}


double topology_charge_continuous(int zi, const VectorField &M, const Field &Ms){
    int dim_x = M.getMesh().nx;
    int dim_y = M.getMesh().ny;

    std::string period; int _i;
    M.getMesh().getPeriodicBC(period, _i);
    bool px = period.find("x") != std::string::npos;
    bool py = period.find("y") != std::string::npos;

    VectorField::const_accessor M_acc(M);
    Field::ro_accessor Ms_acc(Ms);
    double S = 0;   

    for (int i=0;i<dim_x;i++){
        for(int j=0;j<dim_y;j++){
            S += topology_continuous_helper(i, j, zi, px, py, dim_x, dim_y,M_acc,Ms_acc);
        }
    }
    return S/(4*PI);
}

Field topology_charge_density_continuous(const VectorField &M, const Field &Ms)
{
    int dim_x = M.getMesh().nx;
    int dim_y = M.getMesh().ny;
    int dim_z = M.getMesh().nz;
    double delta_x = M.getMesh().dx;
    double delta_y = M.getMesh().dy;
    double delta_z = M.getMesh().dz;

    std::string period; int _i;
    M.getMesh().getPeriodicBC(period, _i);
    bool px = period.find("x") != std::string::npos;
    bool py = period.find("y") != std::string::npos;

    VectorField::const_accessor M_acc(M);
    Field::ro_accessor Ms_acc(Ms);

    Field density(RectangularMesh(dim_x,dim_y,dim_z,delta_x,delta_y,delta_z));
    Field::rw_accessor density_acc(density);

    Vector3d Mij,Mx,My;
    for(int k=0;k<dim_z;k++){
        for (int i=0;i<dim_x;i++){
            for(int j=0;j<dim_y;j++){
                density_acc.at(i,j,k) = topology_continuous_helper(i, j, k, px, py, dim_x, dim_y,M_acc,Ms_acc)/(4*PI*delta_x*delta_y);
            }
        }
    }
    return density;
}

inline double sigmaA(const Vector3d s1,const Vector3d s2,const Vector3d s3){
    return atan2(dot(s1,cross(s2,s3)),1+dot(s1,s2)+dot(s2,s3)+dot(s3,s1))*2;
}


inline double helper_luescher_charge(int i0, int j0, int i1, int j1, int zi, VectorField::const_accessor &M_acc, Field::ro_accessor &Ms_acc){
    double S = 0;
    double Msij = Ms_acc.at(i0,j0,zi);
    double Msij1 = Ms_acc.at(i0,j1,zi);
    double Msi1j = Ms_acc.at(i1,j0,zi);
    double Msi1j1 = Ms_acc.at(i1,j1,zi);

    if(Msij != 0.0 && Msi1j != 0.0 && Msi1j1 != 0.0){
        S += sigmaA(M_acc.get(i0,j0,zi)/Msij,M_acc.get(i1,j0,zi)/Msi1j,M_acc.get(i1,j1,zi)/Msi1j1);
    }

    if(Msij != 0.0 && Msi1j1 != 0.0 && Msij1 != 0.0){
        S += sigmaA(M_acc.get(i0,j0,zi)/Msij,M_acc.get(i1,j1,zi)/Msi1j1,M_acc.get(i0,j1,zi)/Msij1);
    }

    return S;
}

double topology_charge_berg_luescher_dual_lattice(int zi, const VectorField &M, const Field &Ms){
    int dim_x = M.getMesh().nx;
    int dim_y = M.getMesh().ny;

    std::string period; int _i;
    M.getMesh().getPeriodicBC(period, _i);
    bool px = period.find("x") != std::string::npos;
    bool py = period.find("y") != std::string::npos;

    VectorField::const_accessor M_acc(M);
    Field::ro_accessor Ms_acc(Ms);
    double S=0;
    Vector3d Mij,Mx,My;
    for (int i=0;i<dim_x-1;i++){
        for(int j=0;j<dim_y-1;j++){
            S += helper_luescher_charge(i,j,i+1,j+1,zi,M_acc,Ms_acc);
        }
    }
    if(px){
        for(int j=0;j<dim_y-1;j++){
            S += helper_luescher_charge(dim_x-1,j,0,j+1,zi,M_acc,Ms_acc);
        }
    }
    
    if(py){
        for (int i=0;i<dim_x-1;i++){
            S += helper_luescher_charge(i,dim_y-1,i+1,0,zi,M_acc,Ms_acc);
        }
    }

    if(px && py){
        S += helper_luescher_charge(dim_x-1,dim_y-1,0,0,zi,M_acc,Ms_acc);
    }

    return S/(4.0*PI);
}

Field topology_charge_berg_luescher_density_dual_lattice(const VectorField &M, const Field &Ms)
{
    int dim_x = M.getMesh().nx;
    int dim_y = M.getMesh().ny;
    int dim_z = M.getMesh().nz;
    double delta_x = M.getMesh().dx;
    double delta_y = M.getMesh().dy;
    double delta_z = M.getMesh().dz;

    std::string period; int _i;
    M.getMesh().getPeriodicBC(period, _i);
    bool px = period.find("x") != std::string::npos;
    bool py = period.find("y") != std::string::npos;

    VectorField::const_accessor M_acc(M);
    Field::ro_accessor Ms_acc(Ms);

    Field density(RectangularMesh(dim_x,dim_y,dim_z,delta_x,delta_y,delta_z));
    Field::rw_accessor density_acc(density);

    double dS=0;
    for(int k=0;k<dim_z;k++){
        for (int i=0;i<dim_x;i++){
            for(int j=0;j<dim_y;j++){
                dS = 0;
                if(i < dim_x-1 && j < dim_y-1){
                    dS = helper_luescher_charge(i,j,i+1,j+1,k,M_acc,Ms_acc);
                }else if(j==dim_y-1  && i==dim_y-1 && px && py){
                    dS = helper_luescher_charge(dim_x-1,dim_y-1,0,0,k,M_acc,Ms_acc);
                }else if(i==dim_x-1 && px){
                    dS = helper_luescher_charge(dim_x-1,j,0,j+1,k,M_acc,Ms_acc);
                }else if(j==dim_y-1 && py){
                    dS = helper_luescher_charge(i,dim_y-1,i+1,0,k,M_acc,Ms_acc);
                }                
                density_acc.at(i,j,k) = dS/(4*PI*delta_x*delta_y);   
            }
        }
    }
    return density;
}

double topology_charge_berg_luescherhelper(int i, int j, int k, bool px, bool py, int dim_x, int dim_y,VectorField::const_accessor &M_acc,Field::ro_accessor &Ms_acc){
    double Msij = 0;
    double Msi1j = 0;
    double Msij1 = 0;
    double Msim1j = 0;
    double Msijm1 = 0;
    double area = 0;   
    double dS  = 0;

    Vector3d vMsij,vMsi1j,vMsij1,vMsim1j,vMsijm1;
    Vector3d Mij,Mx,My;

    Msij = Ms_acc.at(i,j,k);
    if(Msij != 0.0){
        vMsij = M_acc.get(i,j,k);
    }
    dS = 0;
    area = 0;
    if(i >= 1 && j >= 1 && i < dim_x - 1 && j < dim_y - 1){
        Msij1 = Ms_acc.at(i,j+1,k);
        if(Msij1 != 0.0){
            vMsij1 = M_acc.get(i,j+1,k);
        }
        Msi1j = Ms_acc.at(i+1,j,k);
        if(Msi1j != 0.0){
            vMsi1j = M_acc.get(i+1,j,k);
        }
        Msijm1 = Ms_acc.at(i,j-1,k);
        if(Msijm1 != 0.0){
            vMsijm1 = M_acc.get(i,j-1,k);
        }
        Msim1j = Ms_acc.at(i-1,j,k);
        if(Msim1j != 0.0){
            vMsim1j = M_acc.get(i-1,j,k);
        }
    }else{
        Msij1 = 0.0;
        Msi1j = 0.0;
        Msijm1 = 0.0;
        Msim1j = 0.0;
        //Msij1 = Ms_acc.at(i,j+1,k);
        if(j == dim_y-1){
            if(py){
                Msij1 = Ms_acc.at(i,0,k);
                vMsij1 = M_acc.get(i,0,k);
            }
        }else{
            Msij1 = Ms_acc.at(i,j+1,k);
            vMsij1 = M_acc.get(i,j+1,k);
        }
        //Msi1j = Ms_acc.at(i+1,j,k);
        if(i == dim_x-1){
            if(px){
                Msi1j = Ms_acc.at(0,j,k);
                vMsi1j = M_acc.get(0,j,k);
            }
        }else{
            Msi1j = Ms_acc.at(i+1,j,k);
            vMsi1j = M_acc.get(i+1,j,k);
        }
        //Msijm1 = Ms_acc.at(i,j-1,k);
        if(j == 0){
            if(py){
                Msijm1 = Ms_acc.at(i,dim_y-1,k);
                vMsijm1 = M_acc.get(i,dim_y-1,k);
            }
        }else{
            Msijm1 = Ms_acc.at(i,j-1,k);
            vMsijm1 = M_acc.get(i,j-1,k);
        }
        //Msim1j = Ms_acc.at(i-1,j,k);
        if(i == 0){
            if(px){
                Msim1j = Ms_acc.at(dim_x-1,j,k);
                vMsim1j = M_acc.get(dim_x-1,j,k);
            }
        }else{
            Msim1j = Ms_acc.at(i-1,j,k);
            vMsim1j = M_acc.get(i-1,j,k);
        } 
    }

    if(Msij != 0.0 && Msi1j != 0.0 && Msij1 != 0.0){
        dS += sigmaA(vMsij/Msij,vMsi1j/Msi1j,vMsij1/Msij1);
        area += 0.5;
    }

    if(Msij != 0.0 && Msijm1 != 0.0 && Msi1j != 0.0){
        dS += sigmaA(vMsij/Msij,vMsijm1/Msijm1,vMsi1j/Msi1j);
        area += 0.5;
    }

    if(Msij != 0.0 && Msim1j != 0.0 && Msijm1 != 0.0){
        dS += sigmaA(vMsij/Msij,vMsim1j/Msim1j,vMsijm1/Msijm1);
        area += 0.5;
    }

    if(Msij != 0.0 && Msij1 != 0.0 && Msim1j != 0.0){
        dS += sigmaA(vMsij/Msij,vMsij1/Msij1,vMsim1j/Msim1j);
        area += 0.5;
    }
    
    if(area > 0.1){
        dS = dS/area;
    }
    return dS;
}


double topology_charge_berg_luescher(int zi, const VectorField &M, const Field &Ms)
{
    int dim_x = M.getMesh().nx;
    int dim_y = M.getMesh().ny;

    std::string period; int _i;
    M.getMesh().getPeriodicBC(period, _i);
    bool px = period.find("x") != std::string::npos;
    bool py = period.find("y") != std::string::npos;

    VectorField::const_accessor M_acc(M);
    Field::ro_accessor Ms_acc(Ms);
    double S = 0;
    for (int i=0;i<dim_x;i++){
        for(int j=0;j<dim_y;j++){
            S += topology_charge_berg_luescherhelper(i,j,zi,px,py, dim_x,dim_y,M_acc,Ms_acc);
        }
    }
    
    return S/(4*PI);
}

Field topology_charge_berg_luescher_density(const VectorField &M, const Field &Ms)
{
    int dim_x = M.getMesh().nx;
    int dim_y = M.getMesh().ny;
    int dim_z = M.getMesh().nz;
    double delta_x = M.getMesh().dx;
    double delta_y = M.getMesh().dy;
    double delta_z = M.getMesh().dz;

    std::string period; int _i;
    M.getMesh().getPeriodicBC(period, _i);
    bool px = period.find("x") != std::string::npos;
    bool py = period.find("y") != std::string::npos;

    VectorField::const_accessor M_acc(M);
    Field::ro_accessor Ms_acc(Ms);

    Field density(RectangularMesh(dim_x,dim_y,dim_z,delta_x,delta_y,delta_z));
    Field::rw_accessor density_acc(density);

    for(int k=0;k<dim_z;k++){
        for (int i=0;i<dim_x;i++){
            for(int j=0;j<dim_y;j++){
                density_acc.at(i,j,k) = topology_charge_berg_luescherhelper(i,j,k,px,py,dim_x,dim_y,M_acc,Ms_acc)/(4*PI*delta_x*delta_y);
                
            }
        }
    }
    return density;
}
