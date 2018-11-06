#pragma once

#include "stdafx.h"
#include <string>
#include <vector>
using namespace std;

void incision_model(int& dt, vector<vector<double>>& Q, vector<vector<double>>& bed_z,
	vector<vector<int>>& link, vector<double>& Ds, vector<vector<double>>& bottom_width,
	vector<vector<double>>& top_width, vector<double>& length, vector<vector<double>>& sed_supply,
	double& lambda, string& type, double& nu, vector<double>& tau_c, vector<double>& bank_k,
	int& dt_output, ofstream& zout_file, double& b, ofstream& width_file,
	vector<vector<double>>& height_RB, vector<vector<double>>& height_LB,
	vector<vector<double>>& angle_RB, vector<vector<double>>& angle_LB,
	vector<vector<double>>& toe_height_RB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& toe_angle_RB, vector<vector<double>>& toe_angle_LB,
	vector<double>& cohesion, vector<double>& phi, vector<double>& weight,
	vector<double>& cohesion_toe, vector<double>& phi_toe, vector<double>& weight_toe,
	vector<double>& p_conc, vector<double>& bed_p_conc,
	string& bank_erosion, vector<int>& n_xs, vector<vector<double>>& dx_array,
	boost::multi_array<double, 3>& ps, boost::multi_array<double, 3>& fs,
	vector<vector<double>>& Rc, vector<double>& sinuosity, vector<vector<double>>& n_bends,
	int& meandering, ofstream& dxout_file, vector<vector<double>>& cohesive_z, vector<vector<double>>& bed_tau_c,
	double& omegac_star, ofstream& D50_file, int& MC, vector<double>& n_chnl, vector<double>& n_fp,
	vector<vector<double>>& fp_angle, vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L,
	vector<double>& bank_bedload_prop, vector<double>& bed_bedload_prop,
	string& input_path, ofstream& loading_file, vector<vector<double>>& LB_x, ofstream& geom_file,
	vector<vector<double>>& knick_height, vector<vector<double>>& knick_kd, vector<vector<double>>& knick_x,
	vector<vector<double>>& knick_z, double& transport_factor, double& fluvial_factor, double& cohesive_factor, 
	double& k_factor, int& dt_Q, vector<vector<double>>& bank_armoring, vector<vector<double>>& bank_veg);

void incision_calcs(int& n_nodes, long& i, int& dt, int&day, vector<vector<double>>& slope,
	vector<vector<double>>& bed_z, vector<double>& Ds, vector<double>& Ds_log, string& type, double& nu, double& lambda,
	vector<vector<double>>& bottom_width, vector<vector<double>>& q_array, vector<vector<int>>& link,
	vector<double>& length, vector<double>& sed_supply, int& dt_output,
	double& b, vector<vector<double>>& height_RB, vector<vector<double>>& height_LB,
	vector<vector<double>>& angle_RB, vector<vector<double>>& angle_LB,
	vector<vector<double>>& toe_height_RB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& toe_angle_RB, vector<vector<double>>& toe_angle_LB,
	vector<int>& n_xs,
	vector<vector<double>>& dx_array, double& max_dz, vector<vector<double>>& ds_old,
	int& n_dclass, boost::multi_array<double, 3>& ps, boost::multi_array<double, 3>& fs,
	double& omegac_star, string& limit_fun, const double& a, vector<vector<double>>& cohesive_z,
	vector<vector<double>>& bed_tau_c, vector<vector<double>>& bank_sed, vector<vector<double>>& LB_x,
	vector<vector<double>>& knick_height, vector<vector<double>>& knick_kd, vector<vector<double>>& knick_x,
	vector<vector<double>>& knick_z, vector<double>& bank_bedload_prop, double& bed_sed_out, double& bed_sed_out2,
	vector<vector<double>>& fp_angle, vector<vector<double>>& fp_width_R,
	vector<vector<double>>& fp_width_L, vector<vector<int>>& skip_bank_erosion_L, vector<vector<int>>& skip_bank_erosion_R,
	double& bank_sed_out, double& bed_sed_in, double& dQ_dx_sum, double& cohesive_bed_vol, double& transport_factor,
	double& cohesive_factor, double& k_factor, vector<vector<double>>& cohesive_sed_mass, string& input_path,
	vector<double>& sed_supply_output, vector<double>& bed_bedload_prop, vector<vector<double>>& cohesive_bedload,
	vector<vector<double>>& knick_bedload, vector<vector<double>>& avg_width);

double cohesive_incision(double& omega, double& bed_tau_c, int& dt, double& cohesive_factor, double& k_factor);

void KnickErosion(vector<vector<double>>& knick_height, vector<vector<double>>& knick_kd, vector<vector<double>>& knick_x,
	vector<vector<double>>& knick_z, vector<vector<double>>& bed_z, vector<vector<double>>& Q_chnl, int& n_nodes,
	vector<int>& n_xs, vector<vector<int>>& link, vector<vector<double>>& dx_array, vector<vector<double>>& height_RB,
	vector<vector<double>>& height_LB, vector<vector<double>>& toe_height_RB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& toe_angle_RB, vector<vector<double>>& toe_angle_LB, 
	vector<vector<double>>& avg_width, double& knick_bed_vol, vector<vector<double>>& knick_sed_mass,
	vector<vector<double>>& cohesive_z, vector<vector<double>>& bed_tau_c, vector<double>& bed_bedload_prop,
	vector<vector<double>>& knick_bedload, double& lambda, int& dt, double& bed_sed_out_vol, string& input_path);


