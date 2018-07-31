#ifndef ADD_FUNCTIONS_H
#define ADD_FUNCTIONS_H

#include "stdafx.h"
using namespace std;

double calc_Dx(vector<double>& Ds, vector<double>& ps, int& n_dclass, const double& x,
	vector<double>& Ds_log, string& input_path);

void printOutputs(ofstream& zout_file, vector<vector<double>>& bed_z,
	vector<vector<double>>& width, ofstream& width_file, double& time,
	vector<double>& sinuosity, ofstream& sinuosity_out_file,
	vector<vector<double>>& dx_array, ofstream& dxout_file,
	vector<double>& Ds, vector<double>& Ds_log, boost::multi_array<double, 3>& ps, 
	ofstream& D50_file, vector<int>& n_xs, ofstream& stream_power_file, vector<vector<double>>& slope,
	vector<vector<double>>& q, string& input_path, ofstream& slope_file, ofstream& Rc_file,
	vector<vector<double>>& Rc, vector<vector<double>>& top_width);

void print_bank_outputs(double& time, ofstream& loading_file,
	vector<vector<double>>& bank_sed_mass, vector<vector<double>>& bank_p_mass,
	vector<vector<double>>& cohesive_sed_mass, vector<vector<double>>& knick_sed_mass,
	vector<vector<double>>& bed_p_mass);

void print_knick_outputs(double& time, ofstream& knick_out_file, vector<vector<double>>& knick_x,
	vector<vector<double>>& dx_array);

void print_geometry(vector<vector<double>>& height_LB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& angle_LB, vector<vector<double>>& toe_angle_LB,
	vector<vector<double>>& height_RB, vector<vector<double>>& toe_height_RB,
	vector<vector<double>>& angle_RB, vector<vector<double>>& toe_angle_RB,
	vector<vector<double>>& bottom_width, vector<vector<double>>& fp_angle,
	vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L, vector<vector<double>>& bed_z,
	vector<vector<double>>& LB_x, ofstream& geom_file, double& time);

void print_geometry_final(vector<vector<double>>& height_LB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& angle_LB, vector<vector<double>>& toe_angle_LB,
	vector<vector<double>>& height_RB, vector<vector<double>>& toe_height_RB,
	vector<vector<double>>& angle_RB, vector<vector<double>>& toe_angle_RB,
	vector<vector<double>>& bottom_width, vector<vector<double>>& fp_angle,
	vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L, vector<vector<double>>& bed_z,
	vector<vector<double>>& LB_x, ofstream& geom_file, double& time, vector<int>& n_xs);

double Limiter(double& ri, string& limit_fun);

int find(vector<vector<int>>& link, int value);

double calc_n_bends(double& n_bends, double& reach_length, double& Rc, double& sinuosity);

void create_chnl_geom(double& height_LB, double& toe_height_LB, double& angle_LB,
	double& toe_angle_LB, double& height_RB, double& toe_height_RB, double& angle_RB,
	double& toe_angle_RB, double& bottom_width, double& fp_angle, double& fp_width_R, double& fp_width_L,
	double& bed_z, vector<double>& chnl_x, vector<double>& chnl_y);

double calc_flow_area(vector<double>& chnl_x, vector<double>& chnl_y, double& h, double& area_fp,
	double& area_chnl, double& perim, double& perim_fp, double& perim_chnl, double& ws_width);

double calc_Manning(double& h, vector<double>& chnl_x, vector<double>& chnl_y, double& Q, double& S,
	double& n_chnl, double& n_fp, double& Q_chnl, double& ws_width);

double calc_dz_adj(double& dz, double& toe_height_RB, double& toe_height_LB, double& toe_angle_RB,
	double& toe_angle_LB, double& angle_RB, double& angle_LB, double& width, double& bank_bedload_prop,
	double& fp_angle, double& fp_width_R, double& fp_width_L,
	double& height_RB, double& height_LB, double& bed_z, double& dz_old,
	vector<double>& chnl_x, vector<double>& chnl_y);

void transport_eq_check(const double& omega, const double& q, const double& Ds, double& qs, const string& type);

void reverse_chnl_geom(double& height_LB, double& toe_height_LB, double& angle_LB,
	double& toe_angle_LB, double& height_RB, double& toe_height_RB, double& angle_RB,
	double& toe_angle_RB, double& fp_width_R, double& fp_width_L,
	vector<double>& chnl_x, vector<double>& chnl_y);

long double adjust_chnl(vector<double>& chnl_x, vector<double>& chnl_y, double& dz, const double h_min);


/*void print_summary_file(string& input_path, uint64& time_init, int& dx, int& dt, string& type, string& bank_erosion,
	int& MC);*/

double meandering_eroded_area(double dx, double n_bends, double Rc, double eroded_dist);

void check_files(string& input_path, vector<string>& file_names, string& type);

#endif
