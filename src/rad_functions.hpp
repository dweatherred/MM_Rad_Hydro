#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

//!Initialize Radiaition and Material Properties
double evaluate_opc(const double T, int q){
    if(q==0){ //Shock wave problem
        return  577.35;
    }else if(q==1){ //Marshak Wave
        return 1.0E+6/ pow(T, 3);
    }
}

double evaluate_cv(const double T, double a, int q){
    
    if(q==0){ //Shock wave
        return 1.4467E+12;;
    }else if (q==1){ //Marshak Wave
        return 1.3487E+11;
    }
}


void rad_initialize(std::vector<double> &T0_r, std::vector<double> &E0_r, std::vector<double> &Er_in, double E_total,
                    std::vector<double> &T0_m, std::vector<double> &E0_m, std::vector<double> &Ek,
                    std::vector<double> &Tk, std::vector<double> &opc, std::vector<double> &abs,
                    std::vector<double> &emis, std::vector<double> &rho, std::vector<double> &cv,
                    std::vector<double> &Pr, double a, double c, double dx, double xf, int p){

    int n = abs.size();
    int q;

    if(p==3){
        q=1;
    }else{
        q=0;
    }//q=0 shock wave, q=1 marshak

    for(int j=0; j<n; j++){
        //Temperature
        if(p==0){//Mach 1.2
            if((j+0.5)*dx < xf/2){
                T0_r[j] = 100;
                T0_m[j] = 100;
            }else{
                T0_r[j] = 119.476;
                T0_m[j] = 119.476;
            }
        }else if(p==1){//Mach 3
            if((j+0.5)*dx < xf/2){
                T0_r[j] = 100;
                T0_m[j] = 100;
            }else{
                T0_r[j] = 366.260705;
                T0_m[j] = 366.260705;
            }
        }else if(p==2){//Sod Shock Tube
            T0_r[j] = 0;
            T0_m[j] = 0;
        }else if(p==3){//Marshak
            T0_r[j] = 0.025;
            T0_m[j] = 0.025;
        }

        //Initialize specific heat capacity in each cell
        cv[j] = evaluate_cv(T0_m[j], a, q);

        //Energy
        E0_r[j] = (a * pow(T0_r[j],4)); // erg/cc
        E0_m[j] = rho[j] * cv[j] * T0_m[j]; // erg/cc
        E_total = E0_r[j] + (rho[j] * cv[j] * T0_m[j]);
        Er_in[j] = E0_r[j];

        //Inital opacity in each cell
        opc[j] = evaluate_opc(T0_m[j], q);

        //Abs and Emis
        abs[j] = opc[j] * c * E0_r[j];
        emis[j] = opc[j] * a * c * pow(T0_m[j],4);

        //Ek and Tk for Newton's Method
        Ek[j] = E0_r[j];
        Tk[j] = T0_m[j];

        //Radiation Pressure
        if (p==2){
            Pr[j] = 0;
        }else{
            Pr[j] = E0_r[j]/3.0; 
        }
    } 
}

void moving_shock_rad_init(std::vector<double> &T0_r, std::vector<double> &E0_r, std::vector<double> &Er_in, double E_total,
                           std::vector<double> &T0_m, std::vector<double> &E0_m, std::vector<double> &Ek,
                           std::vector<double> &Tk, std::vector<double> &opc, std::vector<double> &abs,
                           std::vector<double> &emis, std::vector<double> &rho, std::vector<double> &cv,
                           std::vector<double> &Pr, double a, double c, double dx, double xf){
    int n = abs.size();
    int q=0;

    for(int j=0; j<n; j++){
        if((j+0.5)*dx < xf/2){
                T0_r[j] = 100;
                T0_m[j] = 100;
            }else{
                T0_r[j] = 120.5;
                T0_m[j] = 120.5;
        }

        //Initialize specific heat capacity in each cell
        cv[j] = evaluate_cv(T0_m[j], a, q);

        //Energy
        E0_r[j] = (a * pow(T0_r[j],4)); // erg/cc
        E0_m[j] = rho[j] * cv[j] * T0_m[j]; // erg/cc
        E_total = E0_r[j] + (rho[j] * cv[j] * T0_m[j]);
        Er_in[j] = E0_r[j];

        //Inital opacity in each cell
        opc[j] = evaluate_opc(T0_m[j], q);

        //Abs and Emis
        abs[j] = opc[j] * c * E0_r[j];
        emis[j] = opc[j] * a * c * pow(T0_m[j],4);

        //Ek and Tk for Newton's Method
        Ek[j] = E0_r[j];
        Tk[j] = T0_m[j];

        //Radiation Pressure
        Pr[j] = E0_r[j]/3.0; 
    }
}

//! Radiation MMC
void rad_mmc(const std::vector<double> &E0_r, const std::vector<double> &Pr, const std::vector<double> &grad_u,
             const std::vector<double> &lu_Er, const std::vector<double> &vol_old, const std::vector<double> &vol, 
             std::vector<double> &Es_r,  double dt){

    int n_cells = E0_r.size();

    for(int i=0; i<n_cells; i++){
        Es_r[i] = E0_r[i]*vol_old[i] - dt * (lu_Er[i] + Pr[i]*grad_u[i]);
        Es_r[i] = Es_r[i] / vol[i];
    }
}

//!Radiation Solve

// Calculate diffusion coefficients
void setup(std::vector<double> &opc, const std::vector<double> &T0_m, std::vector<double> &D, 
           std::vector<double> &Dp, std::vector<double> &Dm, std::vector<double> &cv,
           const double c, const double a, int p){

    int n=opc.size();
    int q;

    if(p==1){
        q=1;
    }else{
        q=0;
    }//q=0 shock wave, q=1 marshak

    //Setup Loop (Diffusion Coefficients and Opacity), Constant over timestep
    for(int i=0; i<n; i++){
        opc[i] = evaluate_opc(T0_m[i], q); //q=0 shock wave, q=1 marshak
        cv[i] = evaluate_cv(T0_m[i], a, q); //q=0 shock wave, q=1 marshak
        D[i] = c /(3 * opc[i]);
    }

    if(p==0){
        for(int l=0; l<n; l++){
            //+1/2 and -1/2 Diffusion Terms
            if(l==n-1){
                Dp[l] = D[l];
            }else{
                Dp[l] =  2 * D[l] * D[l+1] / (D[l] + D[l+1]);
            }
            if(l==0){
                Dm[l] = D[l];
            }else{
                Dm[l] = 2 * D[l] * D[l-1] / (D[l] + D[l-1]);
            }
        }

    }else{ //For marshak problem
        for(int l=0; l<n; l++){
            //+1/2 and -1/2 Diffusion Terms
            if(l==n-1){
                Dp[l] = D[l];
            }else{
                double Tf = (T0_m[l] + T0_m[l+1]) * 0.5; 
                Dp[l] = c/(3 * evaluate_opc(Tf, q));
            }
            if(l==0){
                Dm[l] = D[l];
            }else{
                double Tf = (T0_m[l] + T0_m[l-1]) * 0.5;
                Dm[l] = c/(3 * evaluate_opc(Tf, q));
            }
        }
    }
}

//Newtons Method
void newtons(const std::vector<double> &Ek, const std::vector<double> &Th, const std::vector<double> &opc, const std::vector<double> &rho, 
             const std::vector<double> &cv, std::vector<double> &Tk, const double dt, const double c, const double a, const int iter){

    int n=Ek.size();
    std::vector<double> dTk(n); std::vector<double> Tk1(n);
    std::vector<double> FT(n); std::vector<double> dTT(n);

    for(int i=0; i<n; i++){
        //newton step
        for(int k=0; k<iter; k++){

            //Calculate FT with "known" Ek and guess "Tk"
            FT[i] = ((rho[i] * cv[i]) * ((Tk[i] - Th[i])/dt)) - (opc[i] * c * Ek[i]) + (opc[i] * a * c * pow(Tk[i],4));

            //Derivative of FT with respect to Tk
            dTT[i] = ((rho[i] * cv[i]) / dt) + (4 * opc[i] * c * a * pow(Tk[i], 3));

            //Solve for dTk
            dTk[i] = -FT[i] / dTT[i];

            //Solve for Tk+1
            Tk1[i] = Tk[i] + dTk[i];

            //Reassign Tk as Tk+1
            Tk[i] = Tk1[i];

            if(fabs(dTk[i]) < 1.0E-5){
                break;
            }
        }
    }
}

//Matrix setup for Implicit 1-D problem with diffusion
void matrix(const std::vector<double> &Tk, const std::vector<double> &Es_r, const std::vector<double> &opc, 
            const std::vector<double> &Dp, const std::vector<double> &Dm, const std::vector<double> &vol,
            std::vector<double> &plus, std::vector<double> &mid, std::vector<double> &minus, std::vector<double> &rs,
            std::vector<double> &x_new, const double dt, const double c, const double a, int p){

    int n=plus.size();

    std::vector<double> dx_i(n,0);
    std::vector<double> dx_i_hf_m(n,0);
    std::vector<double> dx_i_hf_p(n,0);

    //Find dx values

    for(int i=0; i<n; i++){
        dx_i[i] = x_new[i+1] - x_new[i];
    }
    for(int i=0; i<n; i++){
        if(i==0){
            dx_i_hf_m[i] = dx_i[i]/2;
            dx_i_hf_p[i] = (dx_i[i] + dx_i[i+1])/2;
        }else if (i==n-1){
            dx_i_hf_m[i] = (dx_i[i] + dx_i[i-1])/2;
            dx_i_hf_p[i] = dx_i[i]/2;
        }else{
            dx_i_hf_m[i] = (dx_i[i] + dx_i[i-1])/2;
            dx_i_hf_p[i] = (dx_i[i] + dx_i[i+1])/2;
        }
         
    }
    if (p==1){//Marshak wave boundary conditon
        for(int j=0; j<n; j++){
            //upper diagonal
            plus[j] = - (Dp[j] / dx_i[j]) * (1/dx_i_hf_p[j]);
            //main Diagonal
            if(j==0){
                mid[j] = 1/dt + c/(2*dx_i[j]) + ((Dp[j]/dx_i[j]) * (1/dx_i_hf_p[j])) + (opc[j] * c);
            }else{
                mid[j] = 1/dt + ((Dp[j]/dx_i[j]) * (1/dx_i_hf_p[j])) + ((Dm[j]/dx_i[j]) * (1/dx_i_hf_m[j])) + (opc[j] * c);
            }
            //lower diagonal
            minus[j] = -(Dm[j]/dx_i[j]) * (1/dx_i_hf_m[j]);
            //Right Side of Equation
            if(j==0){
                rs[j] = opc[j] * a * c * pow(Tk[j],4) + (c * a * pow(150,4) * 0.5)/dx_i[j] + Es_r[j]/dt;
            }else{
                rs[j] = (opc[j] * a * c * pow(Tk[j],4)) + (Es_r[j]/dt);
            }
        }
    }else{ //Reflected Boundary Condition
        for(int i=0; i<n; i++){
            if(i==0){
                //upper diagonal
                plus[i] = - (Dp[i] / dx_i[i]) * (1/dx_i_hf_p[i]);
                //main Diagonal
                mid[i] = 1/dt + ((Dp[i]/dx_i[i]) * (1/dx_i_hf_p[i])) + (opc[i] * c);
                //lower diagonal
                minus[i] = 0;
                //Right Side of Equation
                rs[i] = (opc[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);

            }else if (i==n-1){
                //upper diagonal
                plus[i] = 0; 
                //main Diagonal
                mid[i] = 1/dt + ((Dm[i]/dx_i[i]) * (1/dx_i_hf_m[i])) + (opc[i] * c);
                //lower diagonal
                minus[i] = - ((Dm[i]/dx_i[i]) * (1/dx_i_hf_m[i]));
                //Right Side of Equation
                rs[i] = (opc[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);

            }else{
                //upper diagonal
                plus[i] = -(Dp[i] / dx_i[i]) * (1/dx_i_hf_p[i]); 
                //main Diagonal
                mid[i] = 1/dt + ((Dp[i]/dx_i[i]) * (1/dx_i_hf_p[i])) + ((Dm[i]/dx_i[i]) * (1/dx_i_hf_m[i])) + (opc[i] * c);
                //lower diagonal
                minus[i] = -((Dm[i]/dx_i[i]) * (1/dx_i_hf_m[i]));
                //Right Side of Equation
                rs[i] = (opc[i] * a * c * pow(Tk[i],4)) + (Es_r[i]/dt);
            }
        }
    }
}

void residual(const std::vector<double> &plus, const std::vector<double> &mid, const std::vector<double> &minus,
              const std::vector<double> &Er, const std::vector<double> &rs, std::vector<double> &res0)
{

    int n=plus.size();
    for(int j=0; j<n; j++){
                
      if( j==0)
      {
        double res(0);
        res = mid[j]*Er[j] + plus[j]*Er[j+1] - rs[j];
        res0[j] = -res;        
      }
      else if(j==n-1)
      {
        double res(0);
        res = minus[j]*Er[j-1]+ mid[j]*Er[j]  - rs[j];
        res0[j] = -res;
      }
      else
      {
        double res(0);
        res = minus[j]*Er[j-1]+ mid[j]*Er[j] + plus[j]*Er[j+1] - rs[j];
        res0[j] = -res;
      }
    }
    
}

//Thomas Algorithm
void thomas_alg(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, 
                std::vector<double> &r, std::vector<double> &sol){
    
    int nx = a.size();
    std::vector<double> tmp(nx,0);
    std::vector<double> tmp2(nx,0);

    tmp[0] = c[0]/b[0];
    for(int i=1;i<nx-1;++i){
        tmp[i] = c[i]/(b[i]-a[i]*tmp[i-1]);
    }

    tmp2[0] = r[0]/b[0];
    for(int i=1;i<nx;++i)
    {
        tmp2[i] = (r[i]-a[i]*tmp2[i-1])/(b[i]-a[i]*tmp[i-1]);

    }
    
    sol[nx-1] = tmp2[nx-1];
    for(int i=nx-2;i>=0;--i){
        sol[i] = tmp2[i]-tmp[i]*sol[i+1];

    }

}

//Energy deposition step
void e_dep(const std::vector<double> &cv, const std::vector<double> &Th, const std::vector<double> &Tk, std::vector<double> &e){

    int n=e.size();

    for(int i=0; i<n; i++){
        e[i] = e[i] + (cv[i] * (Tk[i]- Th[i]));
    }
}

//Post Processing
void rad_reassign(const std::vector<double> &Tk, const std::vector<double> &Ek, const std::vector<double> &opc, const std::vector<double> &rho,
              const std::vector<double> &cv, std::vector<double> &E0_r, std::vector<double> &T0_r, std::vector<double> &E0_m, 
              std::vector<double> &T0_m, std::vector<double> &abs, std::vector<double> &emis, std::vector<double> &Pr, 
              std::vector<double> &Pm, std::vector<double> &P, const double a, const double c, int p){

    int n=E0_r.size(); 

    for(int q=0; q<n; q++){
    //Radiation Terms
        E0_r[q] = Ek[q];
        T0_r[q] = pow((E0_r[q]/a), 0.25);
        //Material Terms
        T0_m[q] = Tk[q];
        E0_m[q] = rho[q] * cv[q] * T0_m[q];

        //Absorption and Emission Terms
        abs[q] = opc[q] * c * E0_r[q];
        emis[q] = opc[q] * a * c * pow(T0_m[q], 4); 
        
        //Radiation Pressure
        if (p==2){
            Pr[q] = 0;
        }else{
            Pr[q] = E0_r[q] / 3.0;
        }
        
        //Total Pressure
        P[q] = Pr[q] + Pm[q];
    }
}