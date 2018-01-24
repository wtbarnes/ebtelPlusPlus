/*
helper.h
General purpose includes to be used everywhere
*/

#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "boost/array.hpp"
#include "../rsp_toolkit/source/xmlreader.h"

// Structure to hold all input parameters
struct Parameters {
  /* Total simulation time (in s) */
  double total_time;
  /* Timestep (in s); when using adaptive solver, the initial timestep*/
  double tau;
  /* Maximum allowed timestep (in s) when using adaptive solver */
  double tau_max;
  /* Loop half length (in cm) */
  double loop_length;
  /* Truncation error tolerance for adaptive solver */
  double adaptive_solver_error;
  /* Safety factor on allowed timestep for adaptive solver */
  double adaptive_solver_safety;
  /* Heat flux saturation limit; 1/6 is a typical value */
  double saturation_limit;
  double turbulent_mean_free_path;
  /* Nominal conductive c1 value; 6.0 is recommended*/
  double c1_cond0;
  /* Nominal radiative c1 value; 0.6 is recommended */
  double c1_rad0;
  /* Force collision frequency to 0.9 s^-1 to emulate single-species fluid */
  bool force_single_fluid;
  /* Switch for radiative loss correction to c1 */
  bool use_c1_loss_correction;
  /* Switch for gravitational correction to c1 */
  bool use_c1_grav_correction;
  /* Switch for power-law radiative loss function */
  bool use_power_law_radiative_losses;
  /* Switch for classical Spitzer conductivity */
  bool use_flux_limiting;
  /* Switch for calculating DEM; if True, runtimes will be much longer */
  bool calculate_dem;
  /* Switch for saving extra equation terms */
  bool save_terms;
  /* Switch for using the adaptive solver option */
  bool use_adaptive_solver;
  /* Path to output file */
  std::string output_filename;
  /* XML node holding DEM calculation parameters */
  tinyxml2::XMLElement * dem_options;
  /* Correction to ion equation of state */
  double boltzmann_correction;
  /* Ion mass correction to account for He abundance */
  double ion_mass_correction;
  /* Ratio of helium to hydrogen  */
  double helium_to_hydrogen_ratio;
  /* Number of grid points */
  size_t N;
};

// Structure to hold all results
struct Results {
  /* Time (in s) */
  std::vector<double> time;
  /* Electron temperature (in K) */
  std::vector<double> temperature_e;
  /* Ion temperature (in K) */
  std::vector<double> temperature_i;
  /* Electron pressure (in dyne cm^-2 s^-1) */
  std::vector<double> pressure_e;
  /* Ion pressure (in dyne cm^-2 s^-1) */
  std::vector<double> pressure_i;
  /* Number density (in cm^-3) */
  std::vector<double> density;
  /* Velocity (in cm s^-1) */
  std::vector<double> velocity;
  /* Heating rate (in erg cm^-3 s^-1) */
  std::vector<double> heat;
};

// Structure to hold equation terms
struct Terms {
  /* Electron heat flux (in erg cm^-2 s^-1)*/
  std::vector<double> f_e;
  /* Ion heat flux (in erg cm^-2 s^-1)*/
  std::vector<double> f_i;
  /* Radiative_loss (in erg cm^3 s^-1) */
  std::vector<double> radiative_loss;
  /* c1 coefficient */
  std::vector<double> c1;
};

// Generic type for state vectors and derivatives
typedef boost::array<double, 5> state_type;

#endif
