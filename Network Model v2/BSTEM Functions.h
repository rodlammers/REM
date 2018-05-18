#pragma once

#include "stdafx.h"
#include <vector>

using namespace std;

void set_bank_geom(vector<vector<double>>& bank_geom, double& height, double& angle,
	double& toe_length, double& toe_angle);

double calc_min_angle(double& base_Z, vector<double>& Phi, vector<vector<double>>& bank_geom,
	vector<double>& bank_layer_Z, vector<double>& bank_layer_X);

double calc_max_angle(double& base_Z, double& base_X, vector<vector<double>>& bank_geom);

double set_bank_intersect_X(vector<vector<double>>& bank_geom, double& base_Z);

vector<double> calc_pore_water_force(double& base_Z, double& failure_angle,
	vector<double>& Phi, double& Phib, vector<vector<double>>& bank_geom, double& FP_fail_Z,
	vector<double>& bank_layer_Z, vector<double>& unit_weight);

double set_water_bank_intersect(vector<vector<double>>& bank_geom, double& WS_Z);

vector<double> calc_FP_shear_intersect(vector<vector<double>>& bank_geom, double& base_Z, double& base_X,
	double& failure_angle, int& FP_point);

vector<double> calc_weight(vector<vector<double>>& bank_geom, vector<double>& FP_intersect, double& base_X,
	double& base_Z, vector<double>& unit_weight, vector<double>& bank_layer_Z, vector<double>& bank_layer_X,
	double& failure_angle, int& FP_point);

double polygon_area(vector<vector<double>>& polygon, int& num_vertices);

void calc_water_force(vector<vector<double>>& bank_geom, double& base_X, double& base_Z,
	vector<double>& confining_force, vector<double>& confining_angle,
	vector<double>& bank_layer_Z, vector<double>& bank_layer_X);

int calc_num_layers(vector<vector<double>>& bank_geom, double& base_Z, vector<double>& bank_intersect_Z);

double compute_FoS(vector<vector<double>>& bank_geom, double& base_Z, double& base_X,
	double& failure_angle, vector<double>& unit_weight, vector<double>& cohesion, vector<double>& Phi,
	vector<double>& bank_layer_Z, vector<double>& bank_layer_X);

void calc_min_FoS(vector<vector<double>>& bank_geom, vector<double>& bank_layer_Z,
	vector<double>& bank_layer_X, vector<double>& cohesion, vector<double>& unit_weight,
	vector<double>& Phi, double& best_FoS, double& best_failure_angle,
	double& best_base_Z);

void bracket_brent(vector<vector<double>>& bank_geom, vector<double>& bank_layer_X,
	vector<double>& bank_layer_Z, vector<double>& cohesion, vector<double>& unit_weight,
	vector<double>& Phi, double& best_failure_angle, double& base_Z,
	double& best_FoS);

double sign_fn(double& x, double& y);

void shift_fn(double& A, double& B, double& C, double& D);