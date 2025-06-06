// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <cassert>
#include <cmath>
#include <memory>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <vector>
#include <fstream>
#include <iomanip>

#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "reconst.h"
#include "eos.h"
#include "evolve.h"
#include "advance.h"

using Util::map_2d_idx_to_1d;
using Util::map_1d_idx_to_2d;
using Util::hbarc;
using Eigen::MatrixXd;
using Eigen::EigenSolver;

Advance::Advance(const EOS &eosIn, const InitData &DATA_in,
                 std::shared_ptr<HydroSourceBase> hydro_source_ptr_in) :
    DATA(DATA_in), eos(eosIn),
    diss_helper(eosIn, DATA_in),
    minmod(DATA_in),
    reconst_helper(eos, DATA_in.echo_level),
    transport_coeffs_(eosIn, DATA_in) {

    hydro_source_terms_ptr = hydro_source_ptr_in;
    flag_add_hydro_source = false;
    if (hydro_source_terms_ptr) {
        if (DATA.Initial_profile == 42) {
            if (hydro_source_terms_ptr->get_number_of_sources() > 0) {
                flag_add_hydro_source = true;
            }
        } else {
            flag_add_hydro_source = true;
        }
    }
}

//! this function evolves one Runge-Kutta step in tau
void Advance::AdvanceIt(const double tau,
                        SCGrid &arena_prev, SCGrid &arena_current,
                        SCGrid &arena_future, const int rk_flag) {
    const int grid_neta = arena_current.nEta();
    const int grid_nx   = arena_current.nX();
    const int grid_ny   = arena_current.nY();

    #pragma omp parallel for collapse(3) schedule(guided)
    for (int ieta = 0; ieta < grid_neta; ieta++)
    for (int ix   = 0; ix   < grid_nx;   ix++  )
    for (int iy   = 0; iy   < grid_ny;   iy++  ) {
        double eta_s_local = - DATA.eta_size/2. + ieta*DATA.delta_eta;
        double x_local     = - DATA.x_size  /2. +   ix*DATA.delta_x;
        double y_local     = - DATA.y_size  /2. +   iy*DATA.delta_y;

        FirstRKStepT(tau, x_local, y_local, eta_s_local,
                     arena_current, arena_future, arena_prev,
                     ix, iy, ieta, rk_flag);

        if (DATA.viscosity_flag == 1) {
            U_derivative u_derivative_helper(DATA, eos);
            u_derivative_helper.MakedU(tau, arena_prev, arena_current,
                                       ix, iy, ieta);
            double theta_local = u_derivative_helper.calculate_expansion_rate(
                                            tau, arena_current, ieta, ix, iy);
            DumuVec a_local;
            u_derivative_helper.calculate_Du_supmu(tau, arena_current,
                                                   ieta, ix, iy, a_local);

            VelocityShearVec sigma_local;
            u_derivative_helper.calculate_velocity_shear_tensor(
                    tau, arena_current, ieta, ix, iy, a_local, sigma_local);

            VorticityVec omega_local;
            u_derivative_helper.calculate_kinetic_vorticity_with_spatial_projector(
                    tau, arena_current, ieta, ix, iy, a_local, omega_local);

            DmuMuBoverTVec baryon_diffusion_vector;
            u_derivative_helper.get_DmuMuBoverTVec(baryon_diffusion_vector);

            FirstRKStepW(tau, arena_prev, arena_current, arena_future, rk_flag,
                         theta_local, a_local, sigma_local, omega_local,
                         baryon_diffusion_vector, ieta, ix, iy);
        }
    }
}


/* %%%%%%%%%%%%%%%%%%%%%% First steps begins here %%%%%%%%%%%%%%%%%% */
void Advance::FirstRKStepT(
        const double tau, const double x_local, const double y_local,
        const double eta_s_local,
        SCGrid &arena_current, SCGrid &arena_future, SCGrid &arena_prev,
        const int ix, const int iy, const int ieta, const int rk_flag) {
    // this advances the ideal part
    double tau_rk = tau + rk_flag*(DATA.delta_tau);

    // Solve partial_a T^{a mu} = -partial_a W^{a mu}
    // Update T^{mu nu}
    // MakeDelatQI gets
    //   qi = q0 if rk_flag = 0 or
    //   qi = q0 + k1 if rk_flag = 1
    // rhs[alpha] is what MakeDeltaQI outputs. 
    // It is the spatial derivative part of partial_a T^{a mu}
    // (including geometric terms)
    TJbVec qi = {0};
    MakeDeltaQI(tau_rk, arena_current, ix, iy, ieta, qi, rk_flag);

    TJbVec qi_source = {0.0};

    if (flag_add_hydro_source) {
        EnergyFlowVec j_mu = {0};
        FlowVec u_local = arena_current(ix,iy,ieta).u;

        hydro_source_terms_ptr->get_hydro_energy_source(
                    tau_rk, x_local, y_local, eta_s_local, u_local, j_mu);
        for (int ii = 0; ii < 4; ii++) {
            qi_source[ii] = tau_rk*j_mu[ii];
            if (isnan(qi_source[ii])) {
                music_message << "qi_source is nan. i = " << ii;
                music_message.flush("error");
                exit(0);
            }
        }

        if (DATA.turn_on_rhob == 1) {
            qi_source[4] = (
                tau_rk*hydro_source_terms_ptr->get_hydro_rhob_source(
                            tau_rk, x_local, y_local, eta_s_local, u_local));
        }
    }

    // now MakeWSource returns partial_a W^{a mu}
    // (including geometric terms)
    TJbVec dwmn ={0.0};
    diss_helper.MakeWSource(tau_rk, arena_current, arena_prev, ix, iy, ieta,
                            dwmn);
    for (int alpha = 0; alpha < 5; alpha++) {
        /* dwmn is the only one with the minus sign */
        qi[alpha] -= dwmn[alpha]*(DATA.delta_tau);

        // add energy moemntum and net baryon density source terms
        qi[alpha] += qi_source[alpha]*DATA.delta_tau;

        // set baryon density back to zero if viscous correction made it
        // non-zero remove/modify if rho_b!=0
        // - this is only to remove the viscous correction that
        // can make rho_b negative which we do not want.
        //if (DATA.turn_on_rhob == 0) {
        //    if (alpha == 4 && std::abs(qi[alpha]) > 1e-12)
        //        qi[alpha] = 0.;
        //}

        /* if rk_flag > 0, we now have q0 + k1 + k2. 
         * So add q0 and multiply by 1/2 */
        qi[alpha] += rk_flag*get_TJb(arena_prev(ix,iy,ieta), alpha, 0)*tau;
        qi[alpha] *= 1./(1. + rk_flag);
    }

    double tau_next = tau + DATA.delta_tau;
    auto grid_rk_t = reconst_helper.ReconstIt_shell(
                                tau_next, qi, arena_current(ix, iy, ieta)); 
    UpdateTJbRK(grid_rk_t, arena_future(ix, iy, ieta));
}


void Advance::FirstRKStepW(const double tau, SCGrid &arena_prev,
                           SCGrid &arena_current, SCGrid &arena_future,
                           const int rk_flag, const double theta_local,
                           const DumuVec &a_local,
                           const VelocityShearVec &sigma_local,
                           const VorticityVec &omega_local,
                           const DmuMuBoverTVec &baryon_diffusion_vector,
                           const int ieta, const int ix, const int iy) {

    auto grid_pt_prev = &(arena_prev(ix, iy, ieta));
    auto grid_pt_c = &(arena_current(ix, iy, ieta));
    auto grid_pt_f = &(arena_future(ix, iy, ieta));

    const double tau_now  = tau + rk_flag*DATA.delta_tau;

    // Solve partial_a (u^a W^{mu nu}) = 0
    // Update W^{mu nu}
    // mu = 4 is the baryon current qmu

    // calculate delta uWmunu
    // need to use u[0][mu], remember rk_flag = 0 here
    // with the KT flux
    // solve partial_tau (u^0 W^{kl}) = -partial_i (u^i W^{kl}
    /* Advance uWmunu */
    double tempf, temps;
    if (DATA.turn_on_shear == 1) {
        for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
            double w_rhs = 0.;
            int mu = 0;
            int nu = 0;
            map_1d_idx_to_2d(idx_1d, mu, nu);
            diss_helper.Make_uWRHS(tau_now, arena_current, ix, iy, ieta,
                                   mu, nu, w_rhs, theta_local, a_local);
            tempf = (
                  (1. - rk_flag)*(grid_pt_c->Wmunu[idx_1d]*grid_pt_c->u[0])
                + rk_flag*(grid_pt_prev->Wmunu[idx_1d]*grid_pt_prev->u[0])
            );
            temps = diss_helper.Make_uWSource(
                    tau_now, grid_pt_c, grid_pt_prev, mu, nu, rk_flag,
                    theta_local, a_local, sigma_local, omega_local);
            tempf += temps*(DATA.delta_tau);
            tempf += w_rhs;
            tempf += rk_flag*((grid_pt_c->Wmunu[idx_1d])*(grid_pt_c->u[0]));
            tempf *= 1./(1. + rk_flag);
            grid_pt_f->Wmunu[idx_1d] = tempf/(grid_pt_f->u[0]);
        }
    } else {
        for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
            grid_pt_f->Wmunu[idx_1d] = 0.0;
        }
    }

    if (DATA.turn_on_bulk == 1) {
        double p_rhs;
        diss_helper.Make_uPRHS(tau_now, arena_current, ix, iy, ieta,
                               &p_rhs, theta_local);
        tempf = ((1. - rk_flag)*(grid_pt_c->pi_b*grid_pt_c->u[0])
                 + rk_flag*(grid_pt_prev->pi_b*grid_pt_prev->u[0]));
        temps = diss_helper.Make_uPiSource(
                tau_now, grid_pt_c, grid_pt_prev, rk_flag,
                theta_local, sigma_local);
        tempf += temps*(DATA.delta_tau);
        tempf += p_rhs;
        tempf += rk_flag*((grid_pt_c->pi_b)*(grid_pt_c->u[0]));
        tempf *= 1./(1. + rk_flag);
        grid_pt_f->pi_b = tempf/(grid_pt_f->u[0]);
    } else {
        grid_pt_f->pi_b = 0.0;
    }

    // CShen: add source term for baryon diffusion
    if (DATA.turn_on_diff == 1) {
        int mu = 4;
        for (int idx_1d = 11; idx_1d < 14; idx_1d++) {
            int nu = idx_1d - 10;
            double w_rhs = diss_helper.Make_uqRHS(
                        tau_now, arena_current, ix, iy, ieta, mu, nu);
            tempf = ((1. - rk_flag)*(grid_pt_c->Wmunu[idx_1d]*grid_pt_c->u[0])
                     + rk_flag*(grid_pt_prev->Wmunu[idx_1d]*grid_pt_prev->u[0]));
            temps = diss_helper.Make_uqSource(
                        tau_now, grid_pt_c, grid_pt_prev, nu, rk_flag,
                        theta_local, a_local, sigma_local, omega_local,
                        baryon_diffusion_vector);
            tempf += temps*(DATA.delta_tau);
            tempf += w_rhs;

            tempf += rk_flag*(grid_pt_c->Wmunu[idx_1d]*grid_pt_c->u[0]);
            tempf *= 1./(1. + rk_flag);

            grid_pt_f->Wmunu[idx_1d] = tempf/(grid_pt_f->u[0]);
        }
    } else {
        for (int idx_1d = 10; idx_1d < 14; idx_1d++) {
            grid_pt_f->Wmunu[idx_1d] = 0.0;
        }
    }

    // re-make Wmunu[3][3] so that Wmunu[mu][nu] is traceless
    grid_pt_f->Wmunu[9] = (
        (2.*(  grid_pt_f->u[1]*grid_pt_f->u[2]*grid_pt_f->Wmunu[5]
             + grid_pt_f->u[1]*grid_pt_f->u[3]*grid_pt_f->Wmunu[6]
             + grid_pt_f->u[2]*grid_pt_f->u[3]*grid_pt_f->Wmunu[8])
         - (grid_pt_f->u[0]*grid_pt_f->u[0] - grid_pt_f->u[1]*grid_pt_f->u[1])
           *grid_pt_f->Wmunu[4]
         - (grid_pt_f->u[0]*grid_pt_f->u[0] - grid_pt_f->u[2]*grid_pt_f->u[2])
           *grid_pt_f->Wmunu[7])
        /(grid_pt_f->u[0]*grid_pt_f->u[0] - grid_pt_f->u[3]*grid_pt_f->u[3]));

    // make Wmunu[i][0] using the transversality
    for (int mu = 1; mu < 4; mu++) {
        tempf = 0.0;
        for (int nu = 1; nu < 4; nu++) {
            int idx_1d = map_2d_idx_to_1d(mu, nu);
            tempf += grid_pt_f->Wmunu[idx_1d]*grid_pt_f->u[nu];
        }
        grid_pt_f->Wmunu[mu] = tempf/(grid_pt_f->u[0]);
    }

    // make Wmunu[0][0]
    tempf = 0.0;
    for (int nu = 1; nu < 4; nu++)
        tempf += grid_pt_f->Wmunu[nu]*grid_pt_f->u[nu];
    grid_pt_f->Wmunu[0] = tempf/(grid_pt_f->u[0]);

    solveEigenvaluesWmunu(grid_pt_f);

    // make qmu[0] using transversality
    tempf = 0.0;
    for (int nu = 1; nu < 4; nu++) {
        int idx_1d = map_2d_idx_to_1d(4, nu);
        tempf += grid_pt_f->Wmunu[idx_1d]*grid_pt_f->u[nu];
    }
    grid_pt_f->Wmunu[10] = DATA.turn_on_diff*tempf/(grid_pt_f->u[0]);

    // If the energy density of the fluid element is smaller than 0.01GeV
    // reduce Wmunu using the QuestRevert algorithm
    if (DATA.Initial_profile != 0 && DATA.Initial_profile != 1) {
        QuestRevert(tau, grid_pt_f, ieta, ix, iy);
        if (DATA.causality_method == 1){
            nCausalityConstraints(grid_pt_f, tau);
        }else if (DATA.causality_method == 2){
            sCausalityConstraints(grid_pt_f, tau);
        }
        if (DATA.turn_on_diff == 1) {
            QuestRevert_qmu(tau, grid_pt_f, ieta, ix, iy);
        }
    }
}
//Wmunu = 1d array
void Advance::solveEigenvaluesWmunu(Cell_small *grid_pt) {

    MatrixXd A = MatrixXd::Zero(4,4);

    A(0,0) = - grid_pt->Wmunu[0];
    A(0,1) = grid_pt->Wmunu[1];
    A(0,2) = grid_pt->Wmunu[2];
    A(0,3) = grid_pt->Wmunu[3];

    A(1,0) = - grid_pt->Wmunu[1];
    A(1,1) = grid_pt->Wmunu[4];
    A(1,2) = grid_pt->Wmunu[5];
    A(1,3) = grid_pt->Wmunu[6];

    A(2,0) = - grid_pt->Wmunu[2];
    A(2,1) = grid_pt->Wmunu[5];
    A(2,2) = grid_pt->Wmunu[7];
    A(2,3) = grid_pt->Wmunu[8];

    A(3,0) = - grid_pt->Wmunu[3];
    A(3,1) = grid_pt->Wmunu[6];
    A(3,2) = grid_pt->Wmunu[8];
    A(3,3) = grid_pt->Wmunu[9];

    EigenSolver <MatrixXd> es;
    es.compute(A, false);
    
    double min = es.eigenvalues().real()[0];
    double max = min;

    for (int i = 1; i < 4; i++){
        double temp = es.eigenvalues().real()[i];
        if (temp > max){
            max = temp;
        }else if (temp < min){
            min = temp;
        }
    }
    grid_pt->Lambdas[0] = min;
    grid_pt->Lambdas[1] = - min - max;
    grid_pt->Lambdas[2] = max;
}

//check causality first, if violated, calculate alpha and store them in an array. (alpha = 1 otherwise), pick the min alpha between 0 and 1
void Advance::nCausalityConstraints(Cell_small *grid_pt, double tau){
    double eps  = grid_pt->epsilon;
    double rhob = grid_pt->rhob;
    double cs2 = eos.get_cs2(eps, rhob);
    double P = eos.get_pressure(eps, rhob);
    double transportPart_n13 = 2. * 1./transport_coeffs_.get_shear_relax_time_factor();
    double viscousPart1_n13 = transport_coeffs_.get_lambda_piPi_coeff();
    double viscousPart2_n13 = - 1./2. * transport_coeffs_.get_tau_pipi_coeff();
    double transportPart_n56 = cs2 + 4./3. * 1./transport_coeffs_.get_shear_relax_time_factor() + 1./transport_coeffs_.get_bulk_relax_time_factor() * (1./3. - cs2) * (1./3. - cs2);
    double viscousPart1_n56 = 2./3.*transport_coeffs_.get_lambda_piPi_coeff() + transport_coeffs_.get_delta_PiPi_coeff() + cs2;
    double viscousPart2_n56 = transport_coeffs_.get_delta_pipi_coeff() + 1./3.*transport_coeffs_.get_tau_pipi_coeff() + transport_coeffs_.get_lambda_Pipi_coeff() * (1./3.- cs2) + cs2;


    double n1 = transportPart_n13 + viscousPart1_n13 * grid_pt->pi_b / (eps + P) + viscousPart2_n13 * abs(grid_pt->Lambdas[0]) / (eps + P);
    double n3 = transportPart_n13 + viscousPart1_n13 * grid_pt->pi_b / (eps + P) + viscousPart2_n13 * grid_pt->Lambdas[2] / (eps + P);
    double n5 = transportPart_n56 + viscousPart1_n56 * grid_pt->pi_b / (eps + P) + viscousPart2_n56 * grid_pt->Lambdas[0] / (eps + P);
    double n6 = (1. - transportPart_n56) + (1. - viscousPart1_n56) * grid_pt->pi_b / (eps + P) + (1. - viscousPart2_n56) * grid_pt->Lambdas[2] / (eps + P);
    std::vector<double> nCondition {n1, n3, n5, n6};

    double minAlp = 1;

    for (int i = 0; i < nCondition.size(); i++){
        double alp = 1;
        if (nCondition[i] < 0){
            switch(i){
                case 0: //n1
                    alp = - transportPart_n13/(viscousPart1_n13 * grid_pt->pi_b / (eps + P) + viscousPart2_n13 * abs(grid_pt->Lambdas[0])/ (eps + P));
                    break;
                case 1: //n3
                    alp = - transportPart_n13/(viscousPart1_n13 * grid_pt->pi_b / (eps + P) + viscousPart2_n13 * grid_pt->Lambdas[2] / (eps + P));
                    break;
                case 2: //n5
                    alp = - transportPart_n56/(viscousPart1_n56 * grid_pt->pi_b / (eps + P) + viscousPart2_n56 * grid_pt->Lambdas[0] / (eps + P));
                    break;
                case 3: //n6
                    alp = - (1. - transportPart_n56)/((1. - viscousPart1_n56) * grid_pt->pi_b / (eps + P) + (1. - viscousPart2_n56) * grid_pt->Lambdas[2] / (eps + P));
                    break;
            }
        }
        if (alp > 0 && alp < minAlp){
            minAlp = alp;
        }else if (alp < 0){
            minAlp = 0;
        } 
    }

    grid_pt->pi_b = grid_pt->pi_b * minAlp;

    for (double& pi : grid_pt->Wmunu){
        pi = pi * minAlp;
    }
    for (double& lam : grid_pt->Lambdas){
        lam = lam * minAlp;
    }

    //write reduction factor and energy density into a file
    if (eps > 0.01){
    std::ofstream myFile("necessary_causality_reduction_factor_wtau.dat", std::fstream::out | std::fstream::app);
    myFile << std::scientific << std::setw(18) << std::setprecision(8) 
       << minAlp << "   " << eps << "   " << tau << std::endl;
    myFile.close();
    }
    
}

bool Advance::BinarySearch(double left, double right, double (*func)(double, Cell_small*, void*), void* pt2object, double& result, Cell_small *grid_pt){
    while(right - left > 1e-4){
        if (right < left){
            std::cout << "Upper beta bound is smaller than the lower beta bound" << std::endl;
            return false;
        }  
        if (func(right, grid_pt, pt2object)*func(left, grid_pt, pt2object) > 0){
            return false;
        }    
        result = (left + right) / 2.;
        if (func(result, grid_pt, pt2object) < 0){
            right = result;
        }else{
            left = result;
        }
    }
    return true;
}

double Advance::Suff5(double beta, Cell_small *grid_pt){
    double eps = grid_pt->epsilon;
    double rhob = grid_pt->rhob;
    double cs2 = eos.get_cs2(eps, rhob);
    double P = eos.get_pressure(eps, rhob);
    double L1 = grid_pt->Lambdas[0] / (eps + P);
    double L3 = grid_pt->Lambdas[2] / (eps + P);
    double Pi = grid_pt->pi_b / (eps + P);
    double s_relax = 1./transport_coeffs_.get_shear_relax_time_factor();
    double b_relax = 1./transport_coeffs_.get_bulk_relax_time_factor() * (1./3. - cs2) * (1./3. - cs2);
    double lam_piPi = transport_coeffs_.get_lambda_piPi_coeff();
    double tau_pipi = transport_coeffs_.get_tau_pipi_coeff();
    double del_PiPi = transport_coeffs_.get_delta_PiPi_coeff();
    double del_pipi = transport_coeffs_.get_delta_pipi_coeff();
    double lam_Pipi = transport_coeffs_.get_lambda_Pipi_coeff();

    return 1. - cs2 - 4./3.*s_relax - b_relax - beta*((cs2 - 1. + 2./3.*lam_piPi + del_PiPi)*Pi + (del_pipi + 1./3.*tau_pipi + lam_Pipi + cs2)*L3 + abs(L1)) 
    - beta*beta*(del_pipi - 1./12.*tau_pipi)*(lam_Pipi + cs2 - 1./12.*tau_pipi)*(L3 + abs(L1))*(L3 + abs(L1))/(1. - s_relax + beta*((1. - 1./2.*lam_piPi)*Pi - abs(L1) - 1./2.*tau_pipi*L3));
}
double Advance::Suff7(double beta, Cell_small *grid_pt){
    double eps = grid_pt->epsilon;
    double rhob = grid_pt->rhob;
    double cs2 = eos.get_cs2(eps, rhob);
    double P = eos.get_pressure(eps, rhob);
    double L1 = grid_pt->Lambdas[0] / (eps + P);
    double L3 = grid_pt->Lambdas[2] / (eps + P);
    double Pi = grid_pt->pi_b / (eps + P);
    double s_relax = 1./transport_coeffs_.get_shear_relax_time_factor();
    double b_relax = 1./transport_coeffs_.get_bulk_relax_time_factor() * (1./3. - cs2) * (1./3. - cs2);
    double lam_piPi = transport_coeffs_.get_lambda_piPi_coeff();
    double tau_pipi = transport_coeffs_.get_tau_pipi_coeff();
    double del_PiPi = transport_coeffs_.get_delta_PiPi_coeff();
    double del_pipi = transport_coeffs_.get_delta_pipi_coeff();
    double lam_Pipi = transport_coeffs_.get_lambda_Pipi_coeff();

    return (s_relax + beta*(1./2.*lam_piPi*Pi - 1./2.*tau_pipi*abs(L1)))*(s_relax + beta*(1./2.*lam_piPi*Pi - 1./2.*tau_pipi*abs(L1))) 
    - beta*beta*(del_pipi - 1./12.*tau_pipi)*(lam_Pipi + cs2 - 1./12.*tau_pipi)*(L3 + abs(L1))*(L3 + abs(L1));
}
double Advance::Suff8(double beta, Cell_small *grid_pt){
    double eps = grid_pt->epsilon;
    double rhob = grid_pt->rhob;
    double cs2 = eos.get_cs2(eps, rhob);
    double P = eos.get_pressure(eps, rhob);
    double L1 = grid_pt->Lambdas[0] / (eps + P);
    double L2 = grid_pt->Lambdas[1] / (eps + P);
    double L3 = grid_pt->Lambdas[2] / (eps + P);
    double Pi = grid_pt->pi_b / (eps + P);
    double s_relax = 1./transport_coeffs_.get_shear_relax_time_factor();
    double b_relax = 1./transport_coeffs_.get_bulk_relax_time_factor() * (1./3. - cs2) * (1./3. - cs2);
    double lam_piPi = transport_coeffs_.get_lambda_piPi_coeff();
    double tau_pipi = transport_coeffs_.get_tau_pipi_coeff();
    double del_PiPi = transport_coeffs_.get_delta_PiPi_coeff();
    double del_pipi = transport_coeffs_.get_delta_pipi_coeff();
    double lam_Pipi = transport_coeffs_.get_lambda_Pipi_coeff();

    return 4./3.*s_relax + b_relax + cs2 + beta*((2./3.*lam_piPi + del_PiPi + cs2)*Pi - (del_pipi + 1./3.*tau_pipi - lam_Pipi + cs2)*abs(L1)) 
    - (1. + beta*(Pi + L2))*(1. + beta*(Pi + L3))/3./(1. + beta*(Pi - abs(L1)))/(1. + beta*(Pi - abs(L1)))*(1. + 2.*s_relax + beta*((1. + lam_piPi)*Pi - abs(Pi) + tau_pipi*L3));
}

double Advance::Suff5_Hook(double beta, Cell_small *grid_pt, void* pt2object){
    Advance* mySelf = (Advance*) pt2object;
    return mySelf->Suff5(beta, grid_pt);
}
double Advance::Suff7_Hook(double beta, Cell_small *grid_pt, void* pt2object){
    Advance* mySelf = (Advance*) pt2object;
    return mySelf->Suff7(beta, grid_pt);
}
double Advance::Suff8_Hook(double beta, Cell_small *grid_pt, void* pt2object){
    Advance* mySelf = (Advance*) pt2object;
    return mySelf->Suff8(beta, grid_pt);
}
void Advance::sCausalityConstraints(Cell_small *grid_pt, double tau){
    double eps = grid_pt->epsilon;
    double rhob = grid_pt->rhob;
    double cs2 = eos.get_cs2(eps, rhob);
    double P = eos.get_pressure(eps, rhob);
    double L1 = grid_pt->Lambdas[0] / (eps + P);
    double L3 = grid_pt->Lambdas[2] / (eps + P);
    double Pi = grid_pt->pi_b / (eps + P);
    double s_relax = 1./transport_coeffs_.get_shear_relax_time_factor();
    double b_relax = 1./transport_coeffs_.get_bulk_relax_time_factor() * (1./3. - cs2) * (1./3. - cs2);
    double lam_piPi = transport_coeffs_.get_lambda_piPi_coeff();
    double tau_pipi = transport_coeffs_.get_tau_pipi_coeff();
    double del_PiPi = transport_coeffs_.get_delta_PiPi_coeff();
    double del_pipi = transport_coeffs_.get_delta_pipi_coeff();
    double lam_Pipi = transport_coeffs_.get_lambda_Pipi_coeff();

    double s1 = 1. - s_relax - L1 + (1. - 1./2.*lam_piPi)*Pi - 1./2.*tau_pipi*L3; //>=0
    double s2 = 2.*s_relax + lam_piPi*Pi - tau_pipi*abs(L1); //>=0 
    double s6 = 1./3.*s_relax + b_relax + cs2 + (1./6.*lam_piPi + del_PiPi + cs2)*Pi + (1./6.*tau_pipi - del_pipi + lam_Pipi - cs2)*abs(L1);//>=0

    std::vector<double> sCondition {s1, s2, s6};
    double minBeta = 1;

    for (int i = 0; i < sCondition.size(); i++){
        double beta = 1;
        if (sCondition[i] < 0){
            switch(i){
                case 0: //s1
                    beta = (s_relax - 1.)/(-abs(L1) + (1. - 1./2.*lam_piPi)*Pi - 1./2.*tau_pipi*L3);
                    break;
                case 1: //s2
                    beta = (-2.*s_relax)/(lam_piPi*Pi - tau_pipi*abs(L1));
                    break;
                case 2: //s6
                    beta = -(1./3.*s_relax + b_relax + cs2)/((1./6.*lam_piPi + del_PiPi + cs2)*Pi + (1./6.*tau_pipi - del_pipi + lam_Pipi - cs2)*abs(L1));
                    break;
            }
        }
        if (beta > 0 && beta < minBeta){
            minBeta = beta;
        }else if (beta < 0){
            minBeta = 0;
        } 
    }
    if (Suff5(minBeta, grid_pt) < 0){
        double result = 0;
        Advance objA(eos, DATA, hydro_source_terms_ptr);
        bool status = BinarySearch(0., minBeta, Advance::Suff5_Hook, (void*) &objA, result, grid_pt);
        if (status){
            minBeta = result;
        }else if(cs2 < 0.15){
            minBeta = 0.;
        }else{
            std::cout << "Suff5 Fails Binary Search" << std::endl;
        }
    }
    if (Suff7(minBeta, grid_pt) < 0){
        double result = 0;
        Advance objB(eos, DATA, hydro_source_terms_ptr);
        bool status = BinarySearch(0., minBeta, Advance::Suff7_Hook, (void*) &objB, result, grid_pt);
        if (status){
            minBeta = result;
        }else{
            std::cout << "Suff7 Fails Binary Search" << std::endl;
        }
    }
    if (Suff8(minBeta, grid_pt) < 0){
        double result = 0;
        Advance objC(eos, DATA, hydro_source_terms_ptr);
        bool status = BinarySearch(0., minBeta, Advance::Suff8_Hook, (void*) &objC, result, grid_pt);
        if (status){
            minBeta = result;
        }else{
            std::cout << "Suff8 Fails Binary Search" << std::endl;
        }
    }
    grid_pt->pi_b = grid_pt->pi_b * minBeta;

    for (double& pi : grid_pt->Wmunu){
        pi = pi * minBeta;
    }
    for (double& lam : grid_pt->Lambdas){
        lam = lam * minBeta;
    }

    //write reduction factor and energy density into a file
    if (eps > 0.01){
    std::ofstream myFile("sufficient_causality_reduction_factor_wtau.dat", std::fstream::out | std::fstream::app);
    myFile << std::scientific << std::setw(18) << std::setprecision(8) 
       << minBeta << "   " << eps << "   " << tau << std::endl;
    myFile.close(); 
    }

}


// update results after RK evolution to grid_pt
void Advance::UpdateTJbRK(const ReconstCell &grid_rk, Cell_small &grid_pt) {
    grid_pt.epsilon = grid_rk.e;
    grid_pt.rhob    = grid_rk.rhob;
    grid_pt.u       = grid_rk.u;
}/* UpdateTJbRK */


//! this function reduce the size of shear stress tensor and bulk pressure
//! in the dilute region to stablize numerical simulations
void Advance::QuestRevert(const double tau, Cell_small *grid_pt,
                          const int ieta, const int ix, const int iy) {
    double eps_scale = 0.1;   // 1/fm^4
    double e_local   = grid_pt->epsilon;
    double rhob      = grid_pt->rhob;

    // regulation factor in the default MUSIC
    // double factor = 300.*tanh(grid_pt->epsilon/eps_scale);
    double xi = 0.05;
    double factor = 10.*DATA.quest_revert_strength*(
                        1./(exp(-(e_local - eps_scale)/xi) + 1.)
                            - 1./(exp(eps_scale/xi) + 1.));
    double factor_bulk = factor;

    double pi_00 = grid_pt->Wmunu[0];
    double pi_01 = grid_pt->Wmunu[1];
    double pi_02 = grid_pt->Wmunu[2];
    double pi_03 = grid_pt->Wmunu[3];
    double pi_11 = grid_pt->Wmunu[4];
    double pi_12 = grid_pt->Wmunu[5];
    double pi_13 = grid_pt->Wmunu[6];
    double pi_22 = grid_pt->Wmunu[7];
    double pi_23 = grid_pt->Wmunu[8];
    double pi_33 = grid_pt->Wmunu[9];

    double pisize = (pi_00*pi_00 + pi_11*pi_11 + pi_22*pi_22 + pi_33*pi_33
         - 2.*(pi_01*pi_01 + pi_02*pi_02 + pi_03*pi_03)
         + 2.*(pi_12*pi_12 + pi_13*pi_13 + pi_23*pi_23));

    double pi_local = grid_pt->pi_b;
    double bulksize = 3.*pi_local*pi_local;

    double p_local = eos.get_pressure(e_local, rhob);
    double eq_size = e_local*e_local + 3.*p_local*p_local;

    // In default MUSIC
    double rho_shear = sqrt(pisize/eq_size)/factor;
    double rho_bulk  = sqrt(bulksize/eq_size)/factor_bulk;

    // Reducing the shear stress tensor
    double rho_shear_max = 0.1;
    if (std::isnan(rho_shear)) {
        for (int mu = 0; mu < 10; mu++) {
            grid_pt->Wmunu[mu] = 0.0;
        }
    } else if (rho_shear > rho_shear_max) {
        if (e_local > eps_scale && DATA.echo_level > 5) {
            music_message << "ieta = " << ieta << ", ix = " << ix
                          << ", iy = " << iy
                          << ", energy density = " << e_local*hbarc
                          << " GeV/fm^3, shear |pi/(epsilon+3*P)| = "
                          << rho_shear;
            music_message.flush("warning");
        }
        for (int mu = 0; mu < 10; mu++) {
            grid_pt->Wmunu[mu] = (rho_shear_max/rho_shear)*grid_pt->Wmunu[mu];
        }
    }

    // Reducing bulk viscous pressure
    double rho_bulk_max = 0.1;
    if (rho_bulk > rho_bulk_max) {
        if (e_local > eps_scale && DATA.echo_level > 5) {
            music_message << "ieta = " << ieta << ", ix = " << ix
                          << ", iy = " << iy
                          << ", energy density = " << e_local*hbarc
                          << " GeV/fm^3, bulk |Pi/(epsilon+3*P)| = "
                          << rho_bulk;
            music_message.flush("warning");
        }
        grid_pt->pi_b = (rho_bulk_max/rho_bulk)*grid_pt->pi_b;
    }
}


//! this function reduce the size of net baryon diffusion current
//! in the dilute region to stablize numerical simulations
void Advance::QuestRevert_qmu(const double tau, Cell_small *grid_pt,
                              const int ieta, const int ix, const int iy) {
    double eps_scale = 0.1;   // in 1/fm^4

    double xi = 0.05;
    double factor = 10.*DATA.quest_revert_strength*(
                            1./(exp(-(grid_pt->epsilon - eps_scale)/xi) + 1.)
                            - 1./(exp(eps_scale/xi) + 1.));

    double q_mu_local[4];
    for (int i = 0; i < 4; i++) {
        // copy the value from the grid
        q_mu_local[i] = grid_pt->Wmunu[10+i];
    }

    // calculate the size of q^\mu
    double q_size = 0.0;
    for (int i = 0; i < 4; i++) {
        double gfac = (i == 0 ? -1.0 : 1.0);
        q_size += gfac*q_mu_local[i]*q_mu_local[i];
    }

    // first check the positivity of q^mu q_mu
    // (in the conversion of gmn = diag(-+++))
    if (q_size < 0.0) {
        music_message << "Advance::QuestRevert_qmu: q^mu q_mu = " << q_size
                      << " < 0!";
        music_message.flush("warning");
        music_message << "Reset it to zero!!!!";
        music_message.flush("warning");
        for (int i = 0; i < 4; i++) {
            int idx_1d = map_2d_idx_to_1d(4, i);
            grid_pt->Wmunu[idx_1d] = 0.0;
        }
    }

    // reduce the size of q^mu according to rhoB
    double e_local = grid_pt->epsilon;
    double rhob_local = grid_pt->rhob;
    double rho_q = sqrt(q_size/(rhob_local*rhob_local))/factor;
    double rho_q_max = 0.1;
    if (rho_q > rho_q_max) {
        if (e_local > eps_scale && DATA.echo_level > 5) {
            music_message << "ieta = " << ieta << ", ix = " << ix
                          << ", iy = " << iy
                          << ", energy density = " << e_local*hbarc
                          << "GeV/fm^3"
                          << ", rhob = " << rhob_local << "1/fm^3"
                          << "-- diffusion |q/rhob| = " << rho_q;
            music_message.flush("warning");
        }
        for (int i = 0; i < 4; i++) {
            grid_pt->Wmunu[10+i] = (rho_q_max/rho_q)*q_mu_local[i];
        }
    }
}


//! This function computes the rhs array. It computes the spatial
//! derivatives of T^\mu\nu using the KT algorithm
void Advance::MakeDeltaQI(const double tau, SCGrid &arena_current,
                          const int ix, const int iy, const int ieta,
                          TJbVec &qi, const int rk_flag) {
    const double delta[4]   = {0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta};
    const double tau_fac[4] = {0.0, tau, tau, 1.0};

    for (int alpha = 0; alpha < 5; alpha++) {
        qi[alpha] = get_TJb(arena_current(ix, iy, ieta), alpha, 0)*tau;
    }

    TJbVec qiphL   = {0.};
    TJbVec qiphR   = {0.};
    TJbVec qimhL   = {0.};
    TJbVec qimhR   = {0.};

    TJbVec rhs     = {0.};
    EnergyFlowVec T_eta_m = {0.};
    EnergyFlowVec T_eta_p = {0.};
    Neighbourloop(arena_current, ix, iy, ieta, NLAMBDAS{
        for (int alpha = 0; alpha < 5; alpha++) {
            const double gphL = qi[alpha];
            const double gphR = tau*get_TJb(p1, alpha, 0);
            const double gmhL = tau*get_TJb(m1, alpha, 0);
            const double gmhR = qi[alpha];
            const double fphL =  0.5*minmod.minmod_dx(gphR, qi[alpha], gmhL);
            const double fphR = -0.5*minmod.minmod_dx(
                                tau*get_TJb(p2, alpha, 0), gphR, qi[alpha]);
            const double fmhL =  0.5*minmod.minmod_dx(
                                qi[alpha], gmhL, tau*get_TJb(m2, alpha, 0));
            const double fmhR = -fphL;
            qiphL[alpha] = gphL + fphL;
            qiphR[alpha] = gphR + fphR;
            qimhL[alpha] = gmhL + fmhL;
            qimhR[alpha] = gmhR + fmhR;
        }

        // for each direction, reconstruct half-way cells
        // reconstruct e, rhob, and u[4] for half way cells
        auto grid_phL = reconst_helper.ReconstIt_shell(tau, qiphL, c);
        auto grid_phR = reconst_helper.ReconstIt_shell(tau, qiphR, c);
        auto grid_mhL = reconst_helper.ReconstIt_shell(tau, qimhL, c);
        auto grid_mhR = reconst_helper.ReconstIt_shell(tau, qimhR, c);

        double aiphL = MaxSpeed(tau, direction, grid_phL);
        double aiphR = MaxSpeed(tau, direction, grid_phR);
        double aimhL = MaxSpeed(tau, direction, grid_mhL);
        double aimhR = MaxSpeed(tau, direction, grid_mhR);

        double aiph = std::max(aiphL, aiphR);
        double aimh = std::max(aimhL, aimhR);

        for (int alpha = 0; alpha < 5; alpha++) {
            double FiphL = get_TJb(grid_phL, 0, alpha, direction)*tau_fac[direction];
            double FiphR = get_TJb(grid_phR, 0, alpha, direction)*tau_fac[direction];
            double FimhL = get_TJb(grid_mhL, 0, alpha, direction)*tau_fac[direction];
            double FimhR = get_TJb(grid_mhR, 0, alpha, direction)*tau_fac[direction];

            // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
            //                  - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
            double Fiph = 0.5*((FiphL + FiphR)
                               - aiph*(qiphR[alpha] - qiphL[alpha]));
            double Fimh = 0.5*((FimhL + FimhR)
                               - aimh*(qimhR[alpha] - qimhL[alpha]));
            if (direction == 3 && (alpha == 0 || alpha == 3)) {
                T_eta_m[alpha] = Fimh;
                T_eta_p[alpha] = Fiph;
            } else {
                double DFmmp = (Fimh - Fiph)/delta[direction];
                rhs[alpha] += DFmmp*(DATA.delta_tau);
            }
        }
    });

    // add longitudinal flux with discretized geometric terms
    double cosh_deta = cosh(delta[3]/2.)/std::max(delta[3], Util::small_eps);
    double sinh_deta = sinh(delta[3]/2.)/std::max(delta[3], Util::small_eps);
    sinh_deta = std::max(0.5, sinh_deta);
    if (DATA.boost_invariant) {
        // if the simulation is boost-invariant,
        // we directly use the limiting value at \Delta eta = 0
        // Longitudinal derivatives should be 0, we set cosh_eta = 0 here
        cosh_deta = 0.0;
        sinh_deta = 0.5;
    }
    rhs[0] += ((  (T_eta_m[0] - T_eta_p[0])*cosh_deta
                - (T_eta_m[3] + T_eta_p[3])*sinh_deta)*DATA.delta_tau);
    rhs[3] += ((  (T_eta_m[3] - T_eta_p[3])*cosh_deta
                - (T_eta_m[0] + T_eta_p[0])*sinh_deta)*DATA.delta_tau);

    // geometric terms
    //rhs[0] -= get_TJb(arena_current(ix, iy, ieta), 3, 3)*DATA.delta_tau;
    //rhs[3] -= get_TJb(arena_current(ix, iy, ieta), 3, 0)*DATA.delta_tau;

    for (int i = 0; i < 5; i++) {
        qi[i] += rhs[i];
    }
}

// determine the maximum signal propagation speed at the given direction
double Advance::MaxSpeed(const double tau, const int direc,
                         const ReconstCell &grid_p) {
    double g[] = {1., 1., 1./tau};

    double utau    = grid_p.u[0];
    double utau2   = utau*utau;
    double ux      = std::abs(grid_p.u[direc]);
    double ut2mux2 = utau2 - ux*ux;

    double eps  = grid_p.e;
    double rhob = grid_p.rhob;

    double vs2 = eos.get_cs2(eps, rhob);
    double num_temp_sqrt = (ut2mux2 - (ut2mux2 - 1.)*vs2)*vs2;
    double num;
    if (num_temp_sqrt >= 0)  {
        num = utau*ux*(1. - vs2) + sqrt(num_temp_sqrt);
    } else {
        double dpde = eos.get_dpde(eps, rhob);
        double p = eos.get_pressure(eps, rhob);
        double h = p+eps;
        if (dpde < 0.001) {
            num = (sqrt(-(h*dpde*h*(dpde*(-1.0 + ut2mux2) - ut2mux2)))
                   - h*(-1.0 + dpde)*utau*ux);
        } else {
          fprintf(stderr,"WARNING: in MaxSpeed. \n");
          fprintf(stderr, "Expression under sqrt in num=%lf. \n", num_temp_sqrt);
          fprintf(stderr,"at value e=%lf. \n",eps);
          fprintf(stderr,"at value p=%lf. \n",p);
          fprintf(stderr,"at value h=%lf. \n",h);
          fprintf(stderr,"at value rhob=%lf. \n",rhob);
          fprintf(stderr,"at value utau=%lf. \n", utau);
          fprintf(stderr,"at value uk=%lf. \n", ux);
          fprintf(stderr,"at value vs^2=%lf. \n", vs2);
          fprintf(stderr,"at value dpde=%lf. \n", eos.get_dpde(eps, rhob));
          fprintf(stderr,"at value dpdrhob=%lf. \n", eos.get_dpdrhob(eps, rhob));
          fprintf(stderr, "MaxSpeed: exiting.\n");
          exit(1);
        }
    }
    double den = utau2*(1. - vs2) + vs2;
    double f = num/std::max(den, Util::small_eps);
    // check for problems
    if (f < 0.0) {
        fprintf(stderr, "SpeedMax = %e\n is negative.\n", f);
        fprintf(stderr, "Can't happen.\n");
        exit(0);
    } else if (f <  ux/utau) {
        if (num != 0.0) {
            if (fabs(f-ux/utau)<0.0001) {
                f = ux/utau;
            } else {
                fprintf(stderr, "SpeedMax-v = %lf\n", f-ux/utau);
                fprintf(stderr, "SpeedMax = %e\n is smaller than v = %e.\n", f, ux/utau);
                fprintf(stderr, "Can't happen.\n");
                exit(0);
            }
        }
    } else if (f > 1.0) {
        fprintf(stderr, "SpeedMax = %e\n is bigger than 1.\n", f);
        fprintf(stderr, "Can't happen.\n");
        fprintf(stderr, "SpeedMax = num/den, num = %e, den = %e \n", num, den);
        fprintf(stderr, "cs2 = %e \n", vs2);
        f =1.;
        exit(1);
    }
    f *= g[direc-1];
    return f;
}

double Advance::get_TJb(const ReconstCell &grid_p, const int rk_flag,
                        const int mu, const int nu) {
    assert(mu < 5); assert(mu > -1);
    assert(nu < 4); assert(nu > -1);
    double rhob = grid_p.rhob;
    const double u_nu = grid_p.u[nu];
    if (mu == 4) {
        return rhob*u_nu;
    }
    double e = grid_p.e;
    double gfac = 0.0;
    double u_mu = 0.0;
    if (mu == nu) {
        u_mu = u_nu;
        gfac = 1.0;
        if (mu == 0) {
            gfac = -1.0;
        }
    } else {
        u_mu = grid_p.u[mu];
    }
    const double pressure = eos.get_pressure(e, rhob);
    const double T_munu   = (e + pressure)*u_mu*u_nu + pressure*gfac;
    return(T_munu);
}

double Advance::get_TJb(const Cell_small &grid_p, const int mu, const int nu) {
    assert(mu < 5); assert(mu > -1);
    assert(nu < 4); assert(nu > -1);
    double rhob = grid_p.rhob;
    const double u_nu = grid_p.u[nu];
    if (mu == 4) {
        return rhob*u_nu;
    }
    double e = grid_p.epsilon;
    double gfac = 0.0;
    double u_mu = 0.0;
    if (mu == nu) {
        u_mu = u_nu;
        gfac = 1.0;
        if (mu == 0) {
            gfac = -1.0;
        }
    } else {
        u_mu = grid_p.u[mu];
    }
    const double pressure = eos.get_pressure(e, rhob);
    const double T_munu   = (e + pressure)*u_mu*u_nu + pressure*gfac;
    return(T_munu);
}
