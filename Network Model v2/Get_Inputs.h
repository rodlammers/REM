#pragma once

#include "stdafx.h"

void input_initial(string& input_type, const string& input_path, int& n_nodes,
	vector<double>& length, vector<int>& n_xs,
	vector<vector<double>>& dx_array, int& max_xs, int& dx);

void input_rest(const string& input_path, string& input_type, int& n_nodes, int& max_xs,
	int& n_dclass, vector<vector<int>>& link, vector<double>& length, vector<int>& n_xs,
	vector<vector<double>>& dx_array, vector<vector<double>>& bed_z, vector<vector<double>>& bottom_width,
	vector<vector<double>>& height_RB, vector<vector<double>>& toe_height_RB, vector<vector<double>>& angle_RB,
	vector<vector<double>>& toe_angle_RB, vector<vector<double>>& height_LB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& angle_LB, vector<vector<double>>& toe_angle_LB, vector<double>& tau_c,
	vector<double>& bank_k, vector<double>& cohesion, vector<double>& phi, vector<double>& weight,
	vector<double>& cohesion_toe, vector<double>& phi_toe, vector<double>& weight_toe, vector<double>& p_conc, 
	vector<double>& bed_p_conc, vector<double>& bank_bedload_prop, vector<double>& bed_bedload_prop,
	vector<vector<double>>& fp_angle, vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L,
	vector<vector<double>>& bed_tau_c, vector<vector<double>>& cohesive_z,
	vector<vector<double>>& top_width, boost::multi_array<double, 3>& ps, boost::multi_array<double, 3>& fs,
	vector<vector<double>>& n_bends, vector<double>& sinuosity, vector<vector<double>>& Rc,
	vector<vector<double>>& LB_x, vector<vector<double>>& knick_height, vector<vector<double>>& knick_kd,
	vector<vector<double>>& knick_x, vector<vector<double>>& knick_z, vector<vector<double>>& bank_armoring,
	vector<vector<double>>& bank_veg);