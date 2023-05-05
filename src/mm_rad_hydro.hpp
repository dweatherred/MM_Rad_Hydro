#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "rad_functions.hpp"
#include "mm_hydro_functions.hpp"

void mm_eularian_rh(double dt, double ti, double tf, double dx, double xl_bound, double xr_bound, int n_cells){

    //Initialize timing
    const int n_time = tf/dt;
    double time = ti;
    double dt_half = dt/2;

    //Number of nodes
    int n_nodes = n_cells + 1;

    //value of gamma
    double mat_gamma = 5.0/3.0;
    double rad_gamma = 1.4;

    //Constants
    const int a = 137;
    double c = 3E+10; // cm/s
    const int n_iter = 1000000000;

    //constant to switch inputs and boundary conditions
    int p;

    //Moving Shock Velocity
    double S;

    //for printing to files
    std::ofstream myfile;

    //Flux terms
    std::vector<double> lu_rho(n_cells,0);
    std::vector<double> lu_v(n_cells,0);
    std::vector<double> lu_e(n_cells,0);
    std::vector<double> lu_w(n_cells,0);
    std::vector<double> lu_Er(n_cells);
    std::vector<double> grad_u(n_cells);

    //!Predictor Values
    //Material and Total values
    std::vector<double> m(n_cells,0); //mass in cell
    std::vector<double> v0(n_cells,0); //initial velocity in cell
    std::vector<double> v_pre(n_cells,0); //updated velocity in cell
    std::vector<double> rho0(n_cells,0); //initial density in cell
    std::vector<double> rho_pre(n_cells,0); //mass in cell
    std::vector<double> P_pre(n_cells,0); //Total pressure in cell
    std::vector<double> Pm_pre(n_cells,0); //Material Pressure in cell
    std::vector<double> e0(n_cells,0); //inital energy in cell
    std::vector<double> e_pre(n_cells,0); //updated energy in cell
    std::vector<double> as(n_cells,0); //speed of sound
    std::vector<double> Th(n_cells, 0); //Hydro Pressure

    //!Corrector Values
    //Material Values
    std::vector<double> v(n_cells,0); //velocity in cell
    std::vector<double> rho(n_cells,0); // density in cell
    std::vector<double> P(n_cells,0); //pressure in cell
    std::vector<double> Pm(n_cells,0); //material pressure in cell
    std::vector<double> e(n_cells,0); // energy in cell

    //!Initial Values That dont get edited
    std::vector<double> v_in(n_cells,0); //velocity in cell
    std::vector<double> rho_in(n_cells,0); // density in cell
    std::vector<double> P_in(n_cells,0); // total pressure in cell
    std::vector<double> e_in(n_cells,0); // energy in cell
    std::vector<double> Er_in(n_cells,0); // Radiation energy in cell
    
    //!Cell Volume
    std::vector<double> vol(n_cells,0); //inital energy in cell
    std::vector<double> vol_old(n_cells,0); //inital energy in cell
    
    //!Moving mesh values
    std::vector<double> w(n_nodes,0); //mesh velocity
    std::vector<double> x0(n_nodes,0); //initital mesh location
    std::vector<double> x_new(n_nodes,0); //updated mesh location

    //!Radiation values
    std::vector<double> T0_r(n_cells);
    std::vector<double> T0_m(n_cells);
    std::vector<double> Es_r(n_cells);
    std::vector<double> E0_r(n_cells);
    std::vector<double> E0_m(n_cells);
    std::vector<double> Pr(n_cells,0); //radiation pressure in cell
    double E_total;

    //! Non predictor or corrector values
    std::vector<double> abs(n_cells);
    std::vector<double> emis(n_cells);
    std::vector<double> opc(n_cells);
    std::vector<double> cv(n_cells);
    std::vector<double> D(n_cells);
    std::vector<double> Dp(n_cells);
    std::vector<double> Dm(n_cells);
    std::vector<double> Ek(n_cells);
    std::vector<double> Tk(n_cells);
    std::vector<double> plus(n_cells);
    std::vector<double> mid(n_cells);
    std::vector<double> minus(n_cells); 
    std::vector<double> rs(n_cells);
    std::vector<double> Er(n_cells);
    std::vector<double> res0(n_cells,0);
    std::vector<double> res1(n_cells,0);    

    //!Initialize Values
    //rad_initialize(T0_r, E0_r, Er_in, E_total, T0_m, E0_m, Ek, Tk, opc, abs, emis, rho, cv, Pr, a, c, dx, xr_bound, p=0); //p=0 mach 1.2, p=1 mach 3, p=2 sod shock tube, p=3 marshak
    //mat_initialize(v0, rho0, Pm, e0, v, e, rho, v_in, e_in, rho_in, P_in, as, Pr, P, T0_m, x0, x_new, m, vol, vol_old, w, dx, mat_gamma, xr_bound, dt, p=0); // p=0 mach 1.2 p=1 mach 3, p=2 sod shock tube, p=3 Marshak

    //!For Moving shock
    moving_shock_rad_init(T0_r, E0_r, Er_in, E_total, T0_m, E0_m, Ek, Tk, opc, abs, emis, rho, cv, Pr, a, c, dx, xr_bound);
    moving_shock_mat_init(v0, rho0, Pm, e0, v, e, rho, v_in, e_in, rho_in, P_in, as, x0, x_new, vol, vol_old, Pr, P, w, T0_m, cv, S, dx, mat_gamma, xr_bound, dt);
    for(int i=0; i<n_time; i++){
        
        //!Hydro Step
        reassign(Pm, v, e, rho, P, w, v0, rho0, e0, as, vol, vol_old, x0, x_new, xl_bound, xr_bound, mat_gamma, dt, p=0); //p=1 for moving shock problems, p=0 for all other moving mesh problems

        //Predictor Flux
        flux(v0, rho0, P, e0, Er, x0, v_in, e_in, rho_in, P_in, Er_in, as, w, lu_rho, lu_v, lu_e, lu_Er, grad_u, vol, mat_gamma, dt_half, p=0); //p=1 for moving shock problems, p=0 for all other moving mesh problems

        //Predictor MM Eularian Calcs
        mm_eularian_calcs(v0, rho0, e0, Pr, lu_rho, lu_v, lu_e, grad_u, cv, vol, vol_old, Pm_pre, v_pre, e_pre, rho_pre, P_pre, Th, as, dt_half, mat_gamma);
        
        //Corrector Flux
        flux(v_pre, rho_pre, P_pre, e_pre, Er, x0, v_in, e_in, rho_in, P_in, Er_in, as, w, lu_rho, lu_v, lu_e, lu_Er, grad_u, vol, mat_gamma, dt, p=0); //p=1 for moving shock problems, p=0 for all other moving mesh problems

        //Corrector MM Eularian Calcs
        mm_eularian_calcs(v0, rho0, e0, Pr, lu_rho, lu_v, lu_e, grad_u, cv, vol, vol_old, Pm, v, e, rho, P, Th, as, dt, mat_gamma);

        //Radiation MMC step
        rad_mmc(E0_r, Pr, grad_u, lu_Er, vol_old, vol, Es_r, dt);

        //Radiation Setup
        setup(opc, T0_m, D, Dp, Dm, cv, c, a, p=0); //p=1 marshak, p=0 else

        //Iteration Loop for implicit temp and energy
        for(int k=0; k<n_iter; k++){

            //Newtons method to solve for temp
            newtons(Ek, Th, opc, rho0, cv, Tk, dt, c, a, n_iter);

            //set up matrix
            matrix(Tk, Es_r, opc, Dp, Dm, vol, plus, mid, minus, rs, x_new, dt, c, a, p=2); //p=1 marshak, p=2 reflective

            //calculate residuals          
            residual(plus, mid, minus, Ek, rs, res0);

            //Thomas Algorithm to solve tridiagonal matrix
            thomas_alg(minus, mid, plus, res0, Er);

            //Check to break out of iteration loop
            double delta_max = -1;
            for(int y=0; y<n_cells; y++){
                double delta = fabs(Er[y]) / Ek[y];
                delta_max = std::max(delta_max, delta);
                Ek[y] += Er[y];
            }
            if(delta_max < 1.0E-5){
                break;
            }
        }

        //Energy Deposition Step
        e_dep(cv, Th, Tk, e);

        //Reassign rad values
        rad_reassign(Tk, Ek, opc, rho, cv, E0_r, T0_r, E0_m, T0_m, abs, emis, Pr, Pm, P, a, c, p=0); //p=2 for sod shock

        double time = dt * i;
        std::cout << " time " << time << std::endl;
    }
    
    myfile.open("../results/MM_moving_shock_500.dat");
    for(int i=0; i<n_cells; i++){   
        myfile << (x0[i]+x0[i+1])*0.5 << " " << T0_r[i] << " " << T0_m[i] << " " << rho[i] << "\n";
    } 
    myfile.close();
}