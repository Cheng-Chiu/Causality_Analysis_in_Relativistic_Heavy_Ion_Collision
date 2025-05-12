# Causality_Analysis_in_Relativistic_Heavy_Ion_Collision

Python notebooks and modified **MUSIC** source files for investigating non-linear causality constraints in 3+1D relativistic viscous hydrodynamics.

# Causality_Project_1

Causality_conditions_checking.ipynb analyzes the non-linear causality constraints presented in Phys. Rev. Lett. 126, 222301. 
hydro_evolution.ipynb applies the non-linear causality constraints to the simulation result generated based on MUSIC, a 3+1D relativistic second-order viscous hydrodynamics for heavy ion collisions written in C++.

# Causality_Project_2

We propose a novel approach, reduction factor method, that dynamically reduce second order viscous variables in-flight during simulatiaon to ensure necessary causality constraints are satisfied. Causality_Checking.ipynb verifies the validity of approach, and Measuring_Reduction_Factor.ipynb measures the evolution of reduction factor during the simulation.

# MUSIC

MUSIC is a 3+1D relativistic second-order viscous hydrodynamics for heavy ion collisions written in C++ (https://github.com/MUSIC-fluid/MUSIC). I incoporated casuality constraint models and corresponding checking in the simulation (advance.cpp, dissipataive.cpp, transport_coeffs.cpp).


