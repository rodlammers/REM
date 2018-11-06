#ifndef BANK_EROSION_H
#define BANK_EROSION_H

#include "stdafx.h"
using namespace std;

double BankShear(const double& omega, const double& width, const double& tau_c, double& bank_k, const int& dt,
	const double& Rc, double& fluvial_factor, double& k_factor, double& bank_armoring);

void BankFailure(double& height, double& angle, double& toe_height, double& toe_angle,
	double& cohesion, double& phi, double& weight, double& cohesion_toe, double& phi_toe,
	double& weight_toe, double& top_width, double& bottom_width, double& fp_angle, double& fp_width,
	double& bank_tank, double& bank_armoring, double& bank_veg);

double FluvialErosion(double& omega, double& toe_angle, double& angle, double& toe_height, double& height,
	double& eroded_area, double& bottom_width, double& tau_c, double& bank_k, int& dt, double& Rc,
	double& top_width, double& fp_width, double& fp_angle, double& bank_tank, int& skip_bank_erosion,
	double& fluvial_factor, double& k_factor, double& bank_armoring);

double calc_eroded_area(double& height, double& toe_height, double& angle,
	double& toe_angle, double& E, double& bottom_width,
	double& top_width, double& fp_width, double& fp_angle);

void FailureBlock(double& toe_angle, double& bottom_width, double& toe_height, 
	double& failblock_area, double& bank_tank);

void BankErosion(vector<double>& tau_c, vector<double>& bank_k, vector<vector<double>>& q,
	vector<vector<double>>& bottom_width, vector<vector<double>>& top_width,
	vector<vector<double>>& slope, int& dt, const int & n_nodes,
	const long& i, vector<vector<double>>& height_RB, vector<vector<double>>& height_LB,
	vector<vector<double>>& angle_RB, vector<vector<double>>& angle_LB, vector<vector<double>>& toe_height_RB,
	vector<vector<double>>& toe_height_LB, vector<vector<double>>& toe_angle_RB,
	vector<vector<double>>& toe_angle_LB, vector<double>& cohesion, vector<double>& phi,
	vector<double>& weight, vector<double>& cohesion_toe, vector<double>& phi_toe,
	vector<double>& weight_toe, string& bank_erosion, vector<int>& n_xs,
	vector<vector<double>>& bank_sed, vector<vector<double>>& dx_array,
	vector<vector<double>>& Rc, vector<double>& sinuosity,
	vector<vector<double>>& n_bends, int& meandering, vector<double>& bank_bedload_prop,
	int& day, int& day_old, vector<vector<double>>& bank_sed_mass, vector<vector<double>>& bank_p_mass,
	vector<double>& p_conc, vector<vector<double>>& LB_x, vector<vector<double>>& fp_angle,
	vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L, vector<vector<double>>& bank_tank_RB,
	vector<vector<double>>& bank_tank_LB, vector<vector<int>>& skip_bank_erosion_L, 
	vector<vector<int>>& skip_bank_erosion_R, double& lambda, double& bank_sed_vol,
	double& fluvial_factor, double& k_factor, vector<vector<double>>& meander_erosion, 
	vector<vector<double>>& bank_armoring, vector<vector<double>>& bank_veg);

#endif
