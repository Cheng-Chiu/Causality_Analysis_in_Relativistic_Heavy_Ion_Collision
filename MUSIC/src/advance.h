// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_ADVANCE_H_
#define SRC_ADVANCE_H_

#include <memory>
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "dissipative.h"
#include "minmod.h"
#include "u_derivative.h"
#include "reconst.h"
#include "hydro_source_base.h"
#include "transport_coeffs.h"
#include "pretty_ostream.h"

class Advance {
 private:
    const InitData &DATA;
    const EOS &eos;
    TransportCoeffs transport_coeffs_;

    std::shared_ptr<HydroSourceBase> hydro_source_terms_ptr;

    Diss diss_helper;
    Minmod minmod;
    Reconst reconst_helper;
    pretty_ostream music_message;

    bool flag_add_hydro_source;

 public:
    Advance(const EOS &eosIn, const InitData &DATA_in,
            std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);

    void AdvanceIt(const double tau_init,
                   SCGrid &arena_prev, SCGrid &arena_current,
                   SCGrid &arena_future, const int rk_flag);

    void FirstRKStepT(const double tau, const double x_local,
                      const double y_local, const double eta_s_local,
                      SCGrid &arena_current, SCGrid &arena_future,
                      SCGrid &arena_prev, const int ix, const int iy,
                      const int ieta, const int rk_flag);

    void FirstRKStepW(const double tau_it, SCGrid &arena_prev,
                      SCGrid &arena_current, SCGrid &arena_future,
                      const int rk_flag, const double theta_local,
                      const DumuVec &a_local,
                      const VelocityShearVec &sigma_local,
                      const VorticityVec &omega_local,
                      const DmuMuBoverTVec &baryon_diffusion_vector,
                      const int ieta, const int ix, const int iy);

    void UpdateTJbRK(const ReconstCell &grid_rk, Cell_small &grid_pt);
    void QuestRevert(const double tau, Cell_small *grid_pt,
                     const int ieta, const int ix, const int iy);
    void QuestRevert_qmu(const double tau, Cell_small *grid_pt,
                         const int ieta, const int ix, const int iy);

    void MakeDeltaQI(const double tau, SCGrid &arena_current,
                     const int ix, const int iy, const int ieta, TJbVec &qi,
                     const int rk_flag);
    void solveEigenvaluesWmunu(Cell_small *grid_pt);
    void nCausalityConstraints(Cell_small *grid_pt, double tau);
    void sCausalityConstraints(Cell_small *grid_pt, double tau);
    bool BinarySearch(double left, double right, double (*func)(double, Cell_small*, void*), void* pt2Object, double& result, Cell_small *grid_pt);
    double Suff5(double beta, Cell_small *grid_pt);
    double Suff7(double beta, Cell_small *grid_pt);
    double Suff8(double beta, Cell_small *grid_pt);
    static double Suff5_Hook(double beta, Cell_small *grid_pt, void* pt2object);
    static double Suff7_Hook(double beta, Cell_small *grid_pt, void* pt2object);
    static double Suff8_Hook(double beta, Cell_small *grid_pt, void* pt2object);

    
    double MaxSpeed(const double tau, const int direc,
                    const ReconstCell &grid_p);

    double get_TJb(const ReconstCell &grid_p, const int rk_flag,
                   const int mu, const int nu);
    double get_TJb(const Cell_small &grid_p, const int mu, const int nu);

};

#endif  // SRC_ADVANCE_H_
