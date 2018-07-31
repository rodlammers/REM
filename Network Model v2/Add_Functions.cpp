
#include "stdafx.h"
#include "Add_Functions.h"

//calc_Dx calculates the grain size of which x * 100% of the distribution is finer (e.g. if
//x = 0.9, this function calculates the D90 of the grain size distribution).

//INPUTS:
//Ds: vector of grain sizes (m)
//ps: vector of grain size fractions (dimensionless; 0-1)
//n_dclass: number of grain sizes
//x: cumulative fraction used to calculate Dx (0-1)
//Ds_log: log2 of grain sizes (to prevent re-calculation every time function is called)

//OUTPUTS:
//Dx: grain size of which x * 100% of the distribution is finer (m)

double calc_Dx(vector<double>& Ds, vector<double>& ps, int& n_dclass, const double& x,
	vector<double>& Ds_log, string& input_path) {
	
	vector<double> ps_cdf(n_dclass); //cumulative distribution of ps

	double Dx; //Dx is what we are solving for

	//If only one grain size, then return the value of that one grain size
	if (n_dclass == 1) {
		Dx = Ds[0];
	}
	else {
		//Get CDF of grain size distribution
		double sum = 0;
		for (int i = 0; i < n_dclass; i++) {
			sum += ps[i];
			ps_cdf[i] = sum;
		}

		//Find left and right bounds for interpolation
		int i_left = -1;
		int i_right = -1;

		//If Dx is smaller than smallest grain size, set to smallest grain size
		if (ps_cdf[0] >= x) {
			Dx = Ds[0];
		}
		else {
			for (int i = 0; i < n_dclass; i++) {
				if (ps_cdf[i] < x) {
					i_left = i;
				}
				if (ps_cdf[i] > x && i_right == -1) {
					i_right = i;
					break;
				}
			}

			//Catch if Dx cannot be calculated for the grain size distribution
			if (i_left < 0 || i_right < 0) {
				cout << i_left << " " << i_right << "\n";
				for (int k = 0; k < n_dclass; k++) {
					cout << ps[k] << " ";
				}
				ofstream error_file(input_path + "MODEL ERRORS.txt");
				error_file << "ERROR IN CALCULATING GRAIN SIZE DISTRIBUTION. CHECK INPUTS OR REDUCE TIME STEP.\n";
				exit(1);
			}

			//Do linear interpolation (in log2 scale)
			double phi_x = Ds_log[i_right] - (ps_cdf[i_right] - x) * (Ds_log[i_right] - Ds_log[i_left]) /
				(ps_cdf[i_right] - ps_cdf[i_left]);

			//Convert back to linear scale (m)
			Dx = pow(2, phi_x);
		}
	}

	return(Dx);

}

//printOutputs: prints outputs of interest to text files at a user-specified interval

//INPUTS:
//zout_file: bed elevation output file
//bed_z: vector of bed elevations to print (m)
//width: vector of channel bottom widths to print (m)
//width_file: width output file
//height: vector of bank heights to print (m)
//heightout_file: height output file
//angle: vector of bank angles to print (radians)
//angleout_file: angle output file
//bank_sed: volumetric loading rate of eroded bank sediment to be added to bed material load (m^3/s)
//bank_sed_out_file: bank_sed output file
//time: simulation day
//sinuosity: reach sinuosity
//sinuosity_out_file: sinuosity output file
//dx_array: channel XS spacing (m) (only changes if modeling meandering)
//dxout_file: dx output file
//Ds: grain sizes (m)
//Ds_log: log2 of grain sizes
//ps: grain size proportions, by XS (0-1)
//D50_file: D50 output file
//n_xs: number of XS per reach

//OUTPUTS: none

void printOutputs(ofstream& zout_file, vector<vector<double>>& bed_z, 
	vector<vector<double>>& width, ofstream& width_file, double& time,
	vector<double>& sinuosity, ofstream& sinuosity_out_file,
	vector<vector<double>>& dx_array, ofstream& dxout_file,
	vector<double>& Ds, vector<double>& Ds_log, boost::multi_array<double, 3>& ps, ofstream& D50_file,
	vector<int>& n_xs, ofstream& stream_power_file, vector<vector<double>>& slope,
	vector<vector<double>>& q, string& input_path, ofstream& slope_file, ofstream& Rc_file,
	vector<vector<double>>& Rc, vector<vector<double>>& top_width) {

	//Find number of nodes (reaches) and max number of XS
	int n_nodes = bed_z.size();
	int n_xs_max = bed_z[0].size();

	//Prints bed elevations for each cross section at the current time step
	for (int j = 0; j < n_nodes - 1; j++) {
		zout_file << time << "	";
		width_file << time << "	";
		dxout_file << time << "	";
		sinuosity_out_file << sinuosity[j] << "	";
		D50_file << time << "	";
		stream_power_file << time << "	";
		slope_file << time << "	";
		Rc_file << time << "	";
		for (int k = 0; k < n_xs_max; ++k) {
			zout_file << bed_z[j][k] << "	";
			width_file << (width[j][k] + top_width[j][k]) / 2.0 << "	";
			dxout_file << dx_array[j][k] << "	";
			stream_power_file << 9810 * slope[j][k] * q[j][k] << "	";
			slope_file << slope[j][k] << "	";
			Rc_file << Rc[j][k] << "	";

			//Calculate D50
			double D50 = 0;
			if (k < n_xs[j]) {
				int n_dclass = Ds.size();
				vector<double> ps_val(n_dclass);
				for (int m = 0; m < n_dclass; ++m) {
					ps_val[m] = ps[j][k][m];
				}
				D50 = calc_Dx(Ds, ps_val, n_dclass, 0.5, Ds_log, input_path) * 1000;
			}

			D50_file << D50 << "	";

		}
		zout_file << "\n";
		width_file << "\n";
		dxout_file << "\n";
		D50_file << "\n";
		stream_power_file << "\n";
		slope_file << "\n";
		Rc_file << "\n";
	}
	sinuosity_out_file << "\n";
}

//print_bank_outputs: prints washload loading from bank erosion at daily timestep (only if bank erosion is modeled)

//INPUTS:
//time: simulation day
//bank_sed_mass_file: output washload mass (kg) file
//bank_p_mass_file: output washload phosphorus (or pollutant) mass (kg) file
//bank_sed_mass: washload mass (kg)
//bank_p_mass: washload phosphorus (or other pollutant) mass (kg)

//OUTPUTS:
//none

void print_bank_outputs(double& time, ofstream& loading_file,
	vector<vector<double>>& bank_sed_mass, vector<vector<double>>& bank_p_mass,
	vector<vector<double>>& cohesive_sed_mass, vector<vector<double>>& knick_sed_mass,
	vector<vector<double>>& bed_p_mass) {

	int n_nodes = bank_sed_mass.size();
	int n_xs_max = bank_sed_mass[0].size();

	//Prints sediment and phosphorus loading for each day
	//Outputs printed as the total loading by reach (summed among all XS in reach)
	//Sediment mass output is sum from bank erosion washload and all cohesive bed erosion and knickpoint erosion
	for (int j = 0; j < n_nodes; j++) {
		//bank_sed_mass_file << time << "	";
		//bank_p_mass_file << time << "	";
		//sed_supply_output_file << sed_supply_output[j] / 86400 << "	"; //Convert m^3/m into rate
		double sed_sum = 0;
		for (int k = 0; k < n_xs_max; ++k) {
			sed_sum += bank_sed_mass[j][k] + cohesive_sed_mass[j][k] + knick_sed_mass[j][k];
		}
		loading_file << sed_sum << "	";
	}
	for (int j = 0; j < n_nodes; j++) {
		double p_sum = 0;
		for (int k = 0; k < n_xs_max; ++k) {
			p_sum += bank_p_mass[j][k] + bed_p_mass[j][k];
		}
		loading_file << p_sum << "	";
	}

	loading_file << "\n";
	//sed_supply_output_file << "\n";

}

void print_knick_outputs(double& time, ofstream& knick_out_file, vector<vector<double>>& knick_x,
	vector<vector<double>>& dx_array) {

	int n_nodes = knick_x.size();
	int n_xs_max = knick_x[0].size();

	for (int j = 0; j < n_nodes; j++) {
		knick_out_file << time << "	";
		double cum_dx = 0;
		for (int k = 0; k < n_xs_max; ++k) {
			double knick_position = 0;
			if (knick_x[j][k] != 0) {
				knick_position = knick_x[j][k] + cum_dx;
			}
			knick_out_file << knick_position << "	";
			cum_dx += dx_array[j][k];
		}
		knick_out_file << "\n";
	}
}
//print_geometry: prints XS geometry for the first XS in each reach at the user-specified time step

//INPUTS:
//height_LB: left bank height (m)
//toe_height_LB: left bank toe height (m)
//angle_LB: left bank angle (radians)
//toe_angle_LB: left bank toe angle (radians)
//height_RB: right bank height (m)
//toe_height_RB: right bank toe height (m)
//angle_RB: right bank angle (radians)
//toe_angle_RB: right bank toe angle (radians)
//bottom_width: channel bottom width (m)
//fp_angle: floodplain angle (radians)
//fp_width_R: right floodplain width (m)
//fp_width_L: left floodplain width (m)
//bed_z: channel bed elevation (m)
//LB_x: station of the bottom of the left bank toe (for plotting) (m)
//geom_file: XS geometry output file
//time: day of simulation

//OUTPUTS:
//none

void print_geometry(vector<vector<double>>& height_LB, vector<vector<double>>& toe_height_LB, 
	vector<vector<double>>& angle_LB, vector<vector<double>>& toe_angle_LB, 
	vector<vector<double>>& height_RB, vector<vector<double>>& toe_height_RB, 
	vector<vector<double>>& angle_RB, vector<vector<double>>& toe_angle_RB, 
	vector<vector<double>>& bottom_width, vector<vector<double>>& fp_angle,
	vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L, vector<vector<double>>& bed_z,
	vector<vector<double>>& LB_x, ofstream& geom_file, double& time) {

	//find number of nodes (reaches)
	int n_nodes = height_LB.size();

	//For each reach, calculate channel geometry and print to file
	for (int j = 0; j < n_nodes; ++j) {
		vector<double> chnl_x(10);
		vector<double> chnl_y(10);
		create_chnl_geom(height_LB[j][0], toe_height_LB[j][0], angle_LB[j][0],
			toe_angle_LB[j][0], height_RB[j][0], toe_height_RB[j][0], angle_RB[j][0],
			toe_angle_RB[j][0], bottom_width[j][0], fp_angle[j][0], fp_width_R[j][0], 
			fp_width_L[j][0], bed_z[j][0],
			chnl_x, chnl_y);

		//print time
		geom_file << time << " ";

		//adjust x values and print
		double diff = LB_x[j][0] - chnl_x[4];
		for (int k = 0; k < 10; ++k) {
			chnl_x[k] += diff;
			geom_file << chnl_x[k] << " ";
		}
		for (int k = 0; k < 10; ++k) {
			geom_file << chnl_y[k] << " ";
		}
		geom_file << "\n";
	}

}

//print_geometry_final: prints the geometry of all XS's at the beginning and end of the simulation

//INPUTS:
//see inputs for "print_geometry" function

//OUTPUTS:
//none

void print_geometry_final(vector<vector<double>>& height_LB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& angle_LB, vector<vector<double>>& toe_angle_LB,
	vector<vector<double>>& height_RB, vector<vector<double>>& toe_height_RB,
	vector<vector<double>>& angle_RB, vector<vector<double>>& toe_angle_RB,
	vector<vector<double>>& bottom_width, vector<vector<double>>& fp_angle,
	vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L, vector<vector<double>>& bed_z,
	vector<vector<double>>& LB_x, ofstream& geom_file, double& time, vector<int>& n_xs) {

	//Find number of nodes (reaches)
	int n_nodes = height_LB.size();

	//For each XS, calculate geometry and print to file
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			vector<double> chnl_x(10);
			vector<double> chnl_y(10);
			create_chnl_geom(height_LB[j][k], toe_height_LB[j][k], angle_LB[j][k],
				toe_angle_LB[j][k], height_RB[j][k], toe_height_RB[j][k], angle_RB[j][k],
				toe_angle_RB[j][k], bottom_width[j][k], fp_angle[j][k], fp_width_R[j][k],
				fp_width_L[j][k], bed_z[j][k],
				chnl_x, chnl_y);

			//print time
			geom_file << time << " ";

			//adjust x values and print
			double diff = LB_x[j][k] - chnl_x[4];
			for (int m = 0; m < 10; ++m) {
				chnl_x[m] += diff;
				geom_file << chnl_x[m] << " ";
			}
			for (int m = 0; m < 10; ++m) {
				geom_file << chnl_y[m] << " ";
			}
			geom_file << "\n";
		}
	}

}

//Limiter: computes the limiter value for use in the flux-limiter sediment transport calculation

//INPUTS:
//ri: ratio of spatial differences in sediment transport rates (m^3/s)
//limit_fun: character string specifying which limiter function to use

//OUTPUTS:
//Psi: value of the limiter

double Limiter(double& ri, string& limit_fun) {

	double Psi;

	if (limit_fun == "Superbee") {
		Psi = max({ 0.0, min(2 * ri, 1.0), min(ri, 2.0) });
	}
	else if (limit_fun == "Van Leer") {
		Psi = (ri + abs(ri)) / (1 + ri);
	}

	return(Psi);
}

//find: Finds index in vector

//INPUTS:
//link: vector linking together nodes (reaches)
//value: value to find in link vector

//OUTPUTS:
//i: column of link vector that value was found in
int find(vector<vector<int>>& link, int value) {
	for (int i = 0; i < link[0].size(); i++) {
		if ((link[0][i] == value) || (link[1][i] == value)) {
			return i;
		}
	}
}

//calc_n_bends: Minimization function to find number of meander bends in a reach given 
//sinuosity and radius of curvature

//INPUTS:
//n_bends: number of meander bends in section (solving for this value)
//reach_length: length of reach (m)
//Rc: bend radius of curvature (m)
//sinuosity: reach sinuosity

//OUTPUT:
//err: calculation error (value to be minimized)

double calc_n_bends(double& n_bends, double& reach_length, double& Rc, double& sinuosity) {
	double bend_length = reach_length / n_bends;
	double angle = bend_length / Rc;
	double Lv = 2 * Rc * sin(0.5 * angle) * n_bends;
	double sinuosity_calc = reach_length / Lv;

	double err = (sinuosity_calc - sinuosity) * (sinuosity_calc - sinuosity);

	return err;
}

//create_chnl_geom: Create channel geometry based on width and bank and floodplain geometry

//INPUTS:
//height_LB: left bank height (m)
//toe_height_LB: left bank toe height (m)
//angle_LB: left bank angle (radians)
//toe_angle_LB: left bank toe angle (radians)
//height_RB: right bank height (m)
//toe_height_RB: right bank toe height (m)
//angle_RB: right bank angle (radians)
//toe_angle_RB: right bank toe angle (radians)
//bottom_width: channel bottom width (m)
//fp_angle: floodplain angle (radians)
//fp_width_R: right floodplain width (m)
//fp_width_L: left floodplain width (m)
//bed_z: channel bed elevation (m)
//chnl_x: vector of XS stations (m)
//chnl_y: vector of XS elevations (m)

//OUTPUT:
//none

void create_chnl_geom(double& height_LB, double& toe_height_LB, double& angle_LB,
	double& toe_angle_LB, double& height_RB, double& toe_height_RB, double& angle_RB,
	double& toe_angle_RB, double& bottom_width, double& fp_angle, double& fp_width_R, double& fp_width_L, 
	double& bed_z, vector<double>& chnl_x, vector<double>& chnl_y) {

	chnl_x[2] = fp_width_L;
	chnl_x[3] = chnl_x[2] + (height_LB - toe_height_LB) / tan(angle_LB);
	chnl_x[4] = chnl_x[3] + toe_height_LB / tan(toe_angle_LB);
	chnl_x[5] = chnl_x[4] + bottom_width;
	chnl_x[6] = chnl_x[5] + toe_height_RB / tan(toe_angle_RB);
	chnl_x[7] = chnl_x[6] + (height_RB - toe_height_RB) / tan(angle_RB);
	chnl_x[8] = chnl_x[7] + fp_width_R;
	chnl_x[9] = chnl_x[8];

	chnl_y[4] = bed_z;
	chnl_y[5] = bed_z;
	chnl_y[3] = chnl_y[4] + toe_height_LB;
	chnl_y[2] = chnl_y[4] + height_LB;
	chnl_y[1] = chnl_y[2] + fp_width_L * tan(fp_angle);
	chnl_y[0] = chnl_y[1] + 20;

	chnl_y[6] = chnl_y[5] + toe_height_RB;
	chnl_y[7] = chnl_y[5] + height_RB;
	chnl_y[8] = chnl_y[7] + fp_width_R * tan(fp_angle);
	chnl_y[9] = chnl_y[8] + 20;
	
	return;
}

//calc_flow_area: Calculate area of flow in channel and overbank

//INPUTS:
//chnl_x: vector of XS stations (m)
//chnl_y: vector of XS elevations (m)
//h: flow depth (m)
//area_fp: flow area in floodplain (m^2)
//area_chnl: flow area in channel (m^2)
//perim: total wetted perimeter (m)
//perim_fp: floodplain wetted perimeter (m)
//perim_chnl: channel wetted perimeter (m)
//ws_width: water surface width (m)

//OUTPUT:
//area: total wetted area (m^2)

double calc_flow_area(vector<double>& chnl_x, vector<double>& chnl_y, double& h, double& area_fp,
	double& area_chnl, double& perim, double& perim_fp, double& perim_chnl, double& ws_width) {

	//Find first and last indices 
	int ind_L = 0, ind_R = 0;
	double WSE = chnl_y[4] + h;
	for (int i = 0; i < 10; ++i) {
		if (chnl_y[i] < WSE) {
			ind_L = i - 1;
			break;
		}
	}

	for (int i = 9; i >= 0; --i) {
		if (chnl_y[i] < WSE) {
			ind_R = i;
			break;
		}
	}

	//Calculate area and perimeter
	double area = 0, sub_area = 0;
	double sub_perim = 0;
	perim = 0, perim_fp = 0, perim_chnl = 0;
	area_fp = 0, area_chnl = 0;
	double y_L = 0, x_L = 0, y_R = 0, x_R = 0;

	for (int i = ind_L; i <= ind_R; ++i) {
		if (i == ind_L) {
			y_L = WSE;
			x_L = chnl_x[i] + (y_L - chnl_y[i]) * (chnl_x[i + 1] - chnl_x[i]) / (chnl_y[i + 1] - chnl_y[i]);
			sub_area = ((WSE - y_L) + (WSE - chnl_y[i + 1])) / 2 * (chnl_x[i + 1] - x_L);
			sub_perim = sqrt((x_L - chnl_x[i + 1]) * (x_L - chnl_x[i + 1]) + (y_L - chnl_y[i + 1]) *
				(y_L - chnl_y[i + 1]));
		}
		else if (i == ind_R) {
			y_R = WSE;
			x_R = chnl_x[i] + (y_R - chnl_y[i]) * (chnl_x[i + 1] - chnl_x[i]) / (chnl_y[i + 1] - chnl_y[i]);
			sub_area = ((WSE - y_R) + (WSE - chnl_y[i])) / 2 * (x_R - chnl_x[i]);
			sub_perim = sqrt((x_R - chnl_x[i]) * (x_R - chnl_x[i]) + (y_R - chnl_y[i]) *
				(y_R - chnl_y[i]));
		}
		else {
			sub_area = ((WSE - chnl_y[i]) + (WSE - chnl_y[i + 1])) / 2 * (chnl_x[i + 1] - chnl_x[i]);
			sub_perim = sqrt((chnl_x[i + 1] - chnl_x[i]) * (chnl_x[i + 1] - chnl_x[i]) + (chnl_y[i + 1] - chnl_y[i]) *
				(chnl_y[i + 1] - chnl_y[i]));
		}

		area += sub_area;
		perim += sub_perim;

		if (i < 2 | i > 6){
			area_fp += sub_area;
			perim_fp += sub_perim;
		}
	}

	area_chnl = area - area_fp;
	perim_chnl = perim - perim_fp;
	ws_width = x_R - x_L;

	return area;
}

//calc_Manning: Minimization function to calculate flow depth and discharge in the channel and overbank

//INPUTS:
//h: flow depth (this is the value we are solving for) (m)
//chnl_x: vector of XS stations (m)
//chnl_y: vector of XS elevations (m)
//Q: discharge (m^3/s)
//S: channel slope
//n_chnl: Manning's n in the channel
//n_fp: Manning's n on the floodplain
//Q_chnl: discharge in the channel (m^3/s)
//ws_width: water surface width (m)

//OUTPUT:
//error: discharge calculation error (this is what is being minimized)

double calc_Manning(double& h, vector<double>& chnl_x, vector<double>& chnl_y, double& Q, double& S,
	double& n_chnl, double& n_fp, double& Q_chnl, double& ws_width) {

	Q_chnl = 0;
	double Q_fp = 0;

	//Full channel flow area, including fp and channel
	double area_fp = 0, area_chnl = 0, perim = 0, perim_fp = 0, perim_chnl = 0;
	double area = calc_flow_area(chnl_x, chnl_y, h, area_fp, area_chnl, perim, perim_fp, perim_chnl,
		ws_width);

	//Channel and fp discharge
	if(area_chnl > 0) Q_chnl = 1 / n_chnl * pow(area_chnl, 5.0 / 3.0) * pow(perim_chnl, -2.0 / 3.0) * sqrt(S);
	if(area_fp > 0) Q_fp = 1 / n_fp * pow(area_fp, 5.0 / 3.0) * pow(perim_fp, -2.0 / 3.0) * sqrt(S);

	double Q_calc = Q_chnl + Q_fp;

	double error = (Q - Q_calc) * (Q - Q_calc);

	return error;
}

//calc_dz_adj: Adjusts magnitude of dz when channel is aggrading (dz > 0) to account for filling up non-rectangular
//channel.

//INPUTS:
//dz: "guessed" value of adjusted dz (m)
//toe_height_RB: right bank toe height (m)
//toe_height_LB: left bank toe height (m)
//toe_angle_RB: right bank toe angle (radians)
//toe_angle_LB: left bank toe angle (radians)
//angle_RB: right bank angle (radians)
//angle_LB: left bank angle (radians)
//width: channel bottom width (m)
//fp_angle: floodplain angle (radians)
//fp_width_R: right floodplain width (m)
//fp_width_L: left floodplain width (m)
//height_RB: right bank height (m)
//height_LB: left bank height (m)
//dz_old: un-adjusted dz from Exner equation (m)
//chnl_x: vector of XS stations (m)
//chnl_y: vector of XS elevations (m)

//OUTPUT:
//error: eroded area error (this is what we are minimizing)

double calc_dz_adj(double& dz, double& toe_height_RB, double& toe_height_LB, double& toe_angle_RB,
double& toe_angle_LB, double& angle_RB, double& angle_LB, double& width, double& bank_bedload_prop,
double& fp_angle, double& fp_width_R, double& fp_width_L,
double& height_RB, double& height_LB, double& bed_z, double& dz_old,
vector<double>& chnl_x, vector<double>& chnl_y) {

	long double A_calc = 0;

	
	//This code section is used if width doesn't change as the channel aggrades
	A_calc = dz * width;
	A_calc += min(toe_height_LB, dz) / 2.0 * (toe_height_LB / tan(toe_angle_LB)) +
		min(toe_height_RB, dz) / 2.0 * (toe_height_RB / tan(toe_angle_RB));

	if (dz > toe_height_RB) {
		//Add area from bank area being filled in
		//A_calc += 0.5 * height_RB * toe_height_RB / tan(toe_angle_RB) - (height_RB - toe_height_RB) *
		//	(toe_height_RB + min(dz, height_RB)) /	(2 * tan(toe_angle_RB)) -
		//	(toe_height_RB * toe_height_RB) / tan(toe_angle_RB) +
		//	(toe_height_RB * min(dz, height_RB)) / (2 * tan(toe_angle_RB));
		A_calc += 0.5 * (min(dz, height_RB) - toe_height_RB) * (min(dz, height_RB) - toe_height_RB) / tan(angle_RB) +
			toe_height_RB / tan(toe_angle_RB) * (min(dz, height_RB) - toe_height_RB);
	}
	if (dz > toe_height_LB) {
		//Add area from bank area being filled in
		//A_calc += 0.5 * height_LB * toe_height_LB / tan(toe_angle_LB) - (height_LB - toe_height_LB) *
		//	(toe_height_LB + min(dz, height_LB)) /	(2 * tan(toe_angle_LB)) -
		//	(toe_height_LB * toe_height_LB) / tan(toe_angle_LB) +
		//	(toe_height_LB * min(dz, height_LB)) / (2 * tan(toe_angle_LB));
		A_calc += 0.5 * (min(dz, height_LB) - toe_height_LB) * (min(dz, height_LB) - toe_height_LB) / tan(angle_LB) +
			toe_height_LB / tan(toe_angle_LB) * (min(dz, height_LB) - toe_height_LB);
	}

	//Adjust areas if whole channel is getting filled in
	if (dz > height_RB) {
		//Add area from floodplain
		if (fp_angle != 0) {
			//A_calc += (2 * (toe_height_RB / tan(toe_angle_RB) + (height_RB - toe_height_RB) / tan(angle_RB)) +
			//	fp_width_R) / 2 * (dz - height_RB);
			//A_banks2 += 0.5 * (dz - height_RB) * (dz - height_RB) / tan(fp_angle);
			A_calc += (toe_height_RB / tan(toe_angle_RB) + (height_RB - toe_height_RB) / tan(angle_RB)) *
				(dz - height_RB) + 0.5 * min((dz - height_RB) / tan(fp_angle), fp_width_R) *
				min(dz - height_RB, fp_width_R * tan(fp_angle)) + max(dz - height_RB - fp_width_R * tan(fp_angle), 0.0) *
				fp_width_R;
		}
		else {
			A_calc += (dz - height_RB) * (fp_width_R + toe_height_RB / tan(toe_angle_RB) +
				(height_RB - toe_height_RB) / tan(angle_RB));
		}
	}
	if (dz > height_LB) {
		//Add area from floodplain
		if (fp_angle != 0) {
			//A_calc += (2 * (toe_height_LB / tan(toe_angle_LB) + (height_LB - toe_height_LB) / tan(angle_LB)) +
			//	fp_width_L) / 2 * (dz - height_LB);
			//A_banks2 += 0.5 * (dz - height_LB) * (dz - height_LB) / tan(fp_angle);
			A_calc += (toe_height_LB / tan(toe_angle_LB) + (height_LB - toe_height_LB) / tan(angle_LB)) *
				(dz - height_LB) + 0.5 * min((dz - height_LB) / tan(fp_angle), fp_width_L) *
				min(dz - height_LB, fp_width_L * tan(fp_angle)) + max(dz - height_LB - fp_width_L * tan(fp_angle), 0.0) *
				fp_width_L;
		}
		else {
			A_calc += (dz - height_LB) * (fp_width_L + toe_height_LB / tan(toe_angle_LB) +
				(height_LB - toe_height_LB) / tan(angle_LB));
		}
	}
	

	//This code section is used of channel widens as channel aggrades
	//Use chnl area calculator to find area filled by dz
	/*double area_fp = 0, area_chnl = 0, perim = 0, perim_fp = 0, perim_chnl = 0, ws_width = 0;
	A_calc = calc_flow_area(chnl_x, chnl_y, dz, area_fp, area_chnl, perim, perim_fp, perim_chnl,
		ws_width);*/

	/*	//New code where width never changes (toe_height, height adjusted to min of 0.01)
		//Need to reset chnl_x and chnl_y every time
	double z_bed = 0;
	vector<double> x(10);
	vector<double> y(10);
	create_chnl_geom(height_LB, toe_height_LB, angle_LB,
		toe_angle_LB, height_RB, toe_height_RB, angle_RB,
		toe_angle_RB, width, fp_angle, fp_width_R,
		fp_width_L, z_bed,
		x, y);
	A_calc = adjust_chnl(x, y, dz, 0.01); */

	double error = (A_calc - abs(dz_old) * width) * (A_calc - abs(dz_old) * width);
	return error;
}

void reverse_chnl_geom(double& height_LB, double& toe_height_LB, double& angle_LB,
	double& toe_angle_LB, double& height_RB, double& toe_height_RB, double& angle_RB,
	double& toe_angle_RB, double& fp_width_R, double& fp_width_L,
	vector<double>& chnl_x, vector<double>& chnl_y) {

	fp_width_L = chnl_x[2] - chnl_x[1];
	fp_width_R = chnl_x[8] - chnl_x[7];

	height_LB = chnl_y[2] - chnl_y[4];
	height_RB = chnl_y[7] - chnl_y[5];

	toe_height_LB = chnl_y[3] - chnl_y[4];
	toe_height_RB = chnl_y[6] - chnl_y[5];

	if ((chnl_x[3] - chnl_x[2]) == 0){
		angle_LB = 90 * M_PI / 180.0;
	}
	else {
		angle_LB = atan((height_LB - toe_height_LB) / (chnl_x[3] - chnl_x[2]));
	}
	if ((chnl_x[7] - chnl_x[6]) == 0) {
		angle_RB = 90 * M_PI / 180.0;
	}
	else {
		angle_RB = atan((height_RB - toe_height_RB) / (chnl_x[7] - chnl_x[6]));
	}

	if ((chnl_x[4] - chnl_x[3]) == 0) {
		toe_angle_LB = 90 * M_PI / 180.0;
	}
	else {
		toe_angle_LB = atan(toe_height_LB / (chnl_x[4] - chnl_x[3]));
	}
	if ((chnl_x[6] - chnl_x[5]) == 0) {
		toe_angle_RB = 90 * M_PI / 180.0;
	}
	else {
		toe_angle_RB = atan(toe_height_RB / (chnl_x[6] - chnl_x[5]));
	}

	if (toe_angle_LB != toe_angle_LB | (toe_angle_RB != toe_angle_RB)) {
		cout << "Toe angle error\n";
	}
}

long double adjust_chnl(vector<double>& chnl_x, vector<double>& chnl_y, double& dz, const double h_min) {
	
	//Store chnl_y and chnl_x values in old vectors
	vector<double> chnl_x_old(10);
	vector<double> chnl_y_old(10);
	for (int i = 0; i < 10; ++i) {
		chnl_x_old[i] = chnl_x[i];
		chnl_y_old[i] = chnl_y[i];
	}
	
	//Find first and last indices 
	int ind_L = 0, ind_R = 0;
	long double h = chnl_y[4] + dz;
	for (int i = 0; i < 10; ++i) {
		if (chnl_y[i] < h) {
			ind_L = i - 1;
			break;
		}
	}

	for (int i = 9; i >= 0; --i) {
		if (chnl_y[i] < h) {
			ind_R = i;
			break;
		}
	}

	double toe_angle_LB = atan((chnl_y_old[3] - chnl_y_old[4]) / (chnl_x_old[4] - chnl_x_old[3]));
	double toe_angle_RB = atan((chnl_y_old[6] - chnl_y_old[5]) / (chnl_x_old[6] - chnl_x_old[5]));
	long double area = dz * (chnl_x[5] - chnl_x[4]);
	area += 0.5 * min(chnl_y[3] - chnl_y[4], dz) * (chnl_y[3] - chnl_y[4]) / tan(toe_angle_LB);
	area += 0.5 * min(chnl_y[6] - chnl_y[5], dz) * (chnl_y[6] - chnl_y[5]) / tan(toe_angle_RB);
	//Adjust bottom elevation
	chnl_y[4] += dz;
	chnl_y[5] += dz;

	if (ind_R == 0) {
		cout << "HELP\n";
	}
	if (ind_R != 0) {
		//Adjust if toes filled in
		if (ind_L < 3) {
			chnl_y[3] = chnl_y[4] + h_min;
			if ((chnl_y_old[2] - chnl_y_old[3]) != 0) {
				chnl_x[3] -= (min(chnl_y[3], chnl_y[2]) - chnl_y_old[3]) * (chnl_x_old[3] - chnl_x_old[2]) /
					(chnl_y_old[2] - chnl_y_old[3]);
			}
			area += (chnl_x_old[4] - chnl_x_old[3]) * 0.5 * ((dz - chnl_y_old[3]) * (chnl_x_old[3] - chnl_x_old[2]) /
				(chnl_y_old[2] - chnl_y_old[3])) + 0.5 * (chnl_x_old[4] - chnl_x_old[3] + ((dz - chnl_y_old[3]) * (chnl_x_old[3] - chnl_x_old[2]) /
				(chnl_y_old[2] - chnl_y_old[3]))) * h_min;

		}
		if (ind_R > 5) {
			chnl_y[6] = chnl_y[4] + h_min;
			if ((chnl_y_old[7] - chnl_y_old[6]) != 0) {
				chnl_x[6] += (min(chnl_y[6], chnl_y[7]) - chnl_y_old[6]) * (chnl_x_old[7] - chnl_x_old[6]) /
					(chnl_y_old[7] - chnl_y_old[6]);
			}
			area += (chnl_x_old[6] - chnl_x_old[5]) * 0.5 * ((dz - chnl_y_old[6]) * (chnl_x_old[7] - chnl_x_old[6]) /
				(chnl_y_old[7] - chnl_y_old[6])) + 0.5 * (chnl_x_old[6] - chnl_x_old[5] + ((dz - chnl_y_old[6]) * (chnl_x_old[7] - chnl_x_old[6]) /
				(chnl_y_old[7] - chnl_y_old[6]))) * h_min;
		}

		//Adjust if banks filled in
		if (ind_L < 2 | chnl_y[2] <= chnl_y[3]) {
			chnl_y[2] = chnl_y[3] + h_min;
			double dist = 0;
			if ((chnl_x_old[2] - chnl_x_old[1]) != 0) {
				dist = (chnl_y[2] - chnl_y_old[2]) / ((chnl_y_old[1] - chnl_y_old[2]) / (chnl_x_old[2] - chnl_x_old[1]));
			}
			if (dist > (chnl_x_old[2] - chnl_x_old[1])) {
				chnl_x[2] = chnl_x[1];
				chnl_y[1] = chnl_y[2];
			}
			else if ((chnl_x_old[2] - chnl_x_old[1]) != 0) {
				chnl_x[2] -= (chnl_y[2] - chnl_y_old[2]) / ((chnl_y_old[1] - chnl_y_old[2]) / (chnl_x_old[2] - chnl_x_old[1]));
			}
		}
		if (ind_R > 6 | chnl_y[7] <= chnl_y[6]) {
			chnl_y[7] = chnl_y[6] + h_min;
			double dist = 0;
			if ((chnl_x_old[8] - chnl_x_old[7]) != 0) {
				dist = (chnl_y[7] - chnl_y_old[7]) / ((chnl_y_old[8] - chnl_y_old[7]) / (chnl_x_old[8] - chnl_x_old[7]));
			}
			if (dist > (chnl_x_old[8] - chnl_x_old[7])) {
				chnl_x[7] = chnl_x[8];
				chnl_y[8] = chnl_y[7];
			}
			else if ((chnl_x_old[8] - chnl_x_old[7]) != 0) {
				chnl_x[7] += (chnl_y[7] - chnl_y_old[7]) / ((chnl_y_old[8] - chnl_y_old[7]) / (chnl_x_old[8] - chnl_x_old[7]));
			}
		}
	}

	if (chnl_x[3] != chnl_x[3]) {
		cout << "ERROR\n";
	}
	if (chnl_x[6] != chnl_x[6]) {
		cout << "ERROR\n";
	}
	//if (ind_R == 0) {
	//	cout << "HELP\n";
	//}
	//Make sure highest points are 20m above edge of fp
	chnl_y[0] = chnl_y[1] + 20;
	chnl_y[9] = chnl_y[8] + 20;

	//Calculate area difference between two XS
	/*vector<long double> A_initial(9), A_final(9);
	for (int i = 1; i < 10; ++i) {
		A_initial[i - 1] = (chnl_x_old[i] - chnl_x_old[i - 1]) * (chnl_y_old[i - 1] + chnl_y_old[i]) / 2;
		A_final[i - 1] = (chnl_x[i] - chnl_x[i - 1]) * (chnl_y[i - 1] + chnl_y[i]) / 2;
	}
	long double area_diff = 0;
	long double c = 0;
	for (int i = 0; i < 9; ++i) {
		long double diff = A_final[i] - A_initial[i];
		//area_diff += diff;
		long double y = diff - c;
		long double t = area_diff + y;
		c = (t - area_diff) - y;
		area_diff = t;
	}
	//double area_diff = A_final - A_initial;
	return(area_diff);*/
	return(area);
}

//transport_eq_check: Checks the input and output values from the sediment transport equations to see if they
//are within the ranges of values used to develop the equations. If not, a warning message is printed at the end
//of the simulation.

//INPUTS:
//omega: specific stream power (W/m^2)
//q: unit discharge (m^2/s)
//Ds: grain size (m)
//qs: sediment transport rate (kg/s for bedload and PPM for total load)
//type: type of transport equation ("bedload" or "total")

//OUTPUT:
//none

void transport_eq_check(const double& omega, const double& q, const double& Ds, double& qs, const string& type) {
	//qs in kg/m/s or ppm for bedload and total load equations, respectively
	if (type == "bedload") {
		if (omega > 412) {
			transport_err_names[0][0] = "Stream power too high.";
			transport_err_vals[0][0] = max(transport_err_vals[0][0], omega);
		}
		if (q > 24.9) {
			transport_err_names[0][1] = "Unit discharge too high.";
			transport_err_vals[0][1] = max(transport_err_vals[0][1], q);
		}
		if (Ds < 0.000137) {
			transport_err_names[0][2] = "Grain size too low.";
			transport_err_vals[0][2] = min(transport_err_vals[0][2], Ds);
		}
		else if (Ds > 0.186) {
			transport_err_names[0][3] = "Grain size too high.";
			transport_err_vals[0][3] = max(transport_err_vals[0][3], Ds);
		}
		if (qs > 12) {
			transport_err_names[0][4] = "Transport rate too high.";
			transport_err_vals[0][4] = max(transport_err_vals[0][4], qs);
		}
	}
	else if (type == "total") {
		if (omega > 161) {
			transport_err_names[1][0] = "Stream power too high";
			transport_err_vals[1][0] = max(transport_err_vals[1][0], omega);
		}
		if (q > 40) {
			transport_err_names[1][1] = "Unit discharge too high";
			transport_err_vals[1][1] = max(transport_err_vals[1][1], q);
		}
		if (Ds < 0.000085) {
			transport_err_names[1][2] = "Grain size too low";
			transport_err_vals[1][2] = min(transport_err_vals[1][2], Ds);
		}
		else if (Ds > 0.0015) {
			transport_err_names[1][3] = "Grain size too high";
			transport_err_vals[1][3] = max(transport_err_vals[1][3], Ds);
		}
		if (qs > 47300) {
			transport_err_names[1][4] = "Transport rate too high";
			transport_err_vals[1][4] = max(transport_err_vals[1][4], qs);
		}
	}

}


/*void print_summary_file(string& input_path, uint64& time_init, int& dx, int& dt, string& type, string& bank_erosion,
	int& MC) {

	//maximum values for sediment transport equations (used in warning message outputs)
	//max specific stream power, max unit discharge, min grain size, max grain size, max transport rate
	vector<vector<double>> equation_range{ { 412, 24.9, 0.000137, 0.186, 12 },
	{ 161, 40, 0.000085, 0.0015, 47300 } };

	ofstream out_file(input_path + "Model Output Summary.txt");

	//Summary of model run
	out_file << "You ran the model with following settings:\n";
	out_file << "XS Spacing: " << dx << " m\n" <<
		"Time step: " << dt << " seconds\n" <<
		"Sediment transport equation: " << type << "\n" <<
		"Bank erosion: " << bank_erosion << "\n" <<
		"Number of Monte Carlo simulations: " << MC << "\n\n";

	//Output the time the model took
	out_file << "Model run time:\n";
	out_file << "Time (seconds): " << (GetTimeMs64() - time_init) / 1000.0 << "\n";
	out_file << "Time (minutes): " << (GetTimeMs64() - time_init) / 1000.0 / 60.0 << "\n";
	out_file << "Time (hours): " << (GetTimeMs64() - time_init) / 1000.0 / 60.0 / 60.0 << "\n";

	out_file << "\n###################################################\n";
	cout << "WARNING MESSAGES\n";
	cout << "\n\n" << "Maximum bed elevation change in time step (abs): " << max_dz << "\n\n";

	double print_ID = 0;
	for (int i = 0; i < 5; ++i) {
		if (transport_err_names[0][i] != "") {
			print_ID = 1;
		}
	}

	if (print_ID == 1) {
		out_file << "You have used the bedload equation outside of the range of values for which it was developed.\n";
		for (int i = 0; i < 5; ++i) {
			if (transport_err_names[0][i] != "") {
				out_file << transport_err_names[0][i] << " Calculated: " << transport_err_vals[0][i] << " Threshold: " <<
					equation_range[0][i] << "\n";
			}
		}
	}

	print_ID = 0;
	for (int i = 0; i < 5; ++i) {
		if (transport_err_names[1][i] != "") {
			print_ID = 1;
		}
	}

	if (print_ID == 1) {
		out_file << "You have used the total load equation outside of the range of values for which it was developed.\n";
		for (int i = 0; i < 5; ++i) {
			if (transport_err_names[1][i] != "") {
				out_file << transport_err_names[1][i] << " Calculated: " << transport_err_vals[1][i] << " Threshold: " <<
					equation_range[1][i] << "\n";
			}
		}
	}

	out_file << flush;

	//Open file for user
	ShellExecuteW(NULL, "open", input_path + "Model Output Summary.txt", "", input_path, SW_SHOW);
}*/

double meandering_eroded_area(double dx, double n_bends, double Rc, double eroded_dist) {

	bool angle_adj = false;

	//Calculate angle of circular bend
	double bend_angle = dx / n_bends / Rc;

	//Calculate segment length
	double seg_length = 2.0 * Rc * sin(0.5 * bend_angle);

	//adjust angle if angle greater than 180 deg
	if (bend_angle > M_PI) { angle_adj = true; }

	//Calculate initial area of curve
	double A_initial = 0.5 * Rc * Rc * (bend_angle - sin(bend_angle));

	//Calculate bend height
	double h = Rc - 0.5 * sqrt(4 * Rc * Rc - seg_length * seg_length);
	//Adjust h depending on the angle
	if (angle_adj) {
		h -= eroded_dist;
	}
	else {
		h += eroded_dist;
	}
	//Calculate new radius of curvature
	Rc = (4 * h * h + seg_length * seg_length) / (8 * h);

	//Calculate new bend angle
	bend_angle = 2 * asin(seg_length / (2 * Rc));

	//Adjust angle if original was greater than 180 deg
	if (angle_adj) { bend_angle = 2.0 * M_PI - bend_angle; }

	//Calculate final area of curve
	double A_final = 0.5 * Rc * Rc * (bend_angle - sin(bend_angle));

	double A_diff = A_final - A_initial;

	return(A_diff);
}

void check_files(string& input_path, vector<string>& file_names, string& type) {

	double n_files = file_names.size();

	for (int i = 0; i < n_files; ++i) {
		ifstream input_file(input_path + file_names[i]);

		if (!input_file) {
			cout << "WARNING: " << type << " INPUT FILE IS MISSING: " << file_names[i] << "\n";

			if (type == "REQUIRED") {
				ofstream error_file(input_path + "MODEL ERRORS.txt");
				error_file << "WARNING: " << type << " INPUT FILE IS MISSING: " << file_names[i] << "\n";
				exit(1);
			}
		}
	}

}