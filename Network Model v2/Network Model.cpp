// Network Model.cpp : Defines the entry point for the console application.
//

#include<iostream>
#include<fstream>
#include<math.h>
#include<string>
#include<vector>
#include<array>
#include<stdlib.h>
#include "stdafx.h"
#include<numeric>
using namespace std;

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

//GLOBAL VARIABLES
vector<vector<string>> transport_err_names(2, vector<string>(5));
vector<vector<double>> transport_err_vals(2, vector<double>(5));

/* Remove if already defined */
typedef long long int64; typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
* windows and linux. */
uint64 GetTimeMs64()
{
#ifdef _WIN32
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	* to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	uint64 ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

	return ret;
#else
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	uint64 ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
	ret /= 1000;

	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);http://127.0.0.1:22347/graphics/plot_zoom_png?width=799&height=900

	return ret;
#endif
}

void print_summary_file(string& input_path, uint64& time_init, int& dx, int& dt, string& type, string& bank_erosion,
	int& MC);

int main() {
	double lambda = 0.4; //bed porosity
	double nu = 0.7; //hiding coefficient - Alonso and Langendoen (0-1)
	double b = 0; //hiding coefficient - new approach (0-(-1.5)) (user defined)
	double omegac_star = 0; //dimensionless critical specific stream power (user defined)
	string type; //total or bedload for sediment transport (or both, EC for Eaton and Church, Parker for Parker)
	string bank_erosion; //none, fluvial, failure, or both
	string input_type; //reach: populate bed_z based on reach elevation, length, and dx
					   //profile: read in full bed profiles and dx's

	int dx, dt, dt_Q; //Spatial step (dx) and time step (dt) and time step of flow data (dt_Q)
	int dt_output; //Frequency of output (number of days)
	int meandering; //incorporate meandering (0 = No, 1 = Yes)
	int MC; //run Monte Carlo simulations (0 = No, anythin else = number of iterations)
	int n_dclass = 0; //number of grain size classes
	int max_threads = 0; //maximum number of threads for parallel processing (usually doesn't need to be set)
	
	//Input/output file locations
	//Check if file exists and set path accordingly
	string input_path;
	ifstream infile("Model Inputs.txt");
	bool exists = infile.good();
	if (exists) {
		input_path = "";
	}
	else {
		input_path = "U:/Network Model v2/";
	}
	infile.close();

	std::cout << "Input file path: " << input_path << "\n";
	
	//Get inputs from file
	ifstream Infile;
	Infile.open(input_path + "Model Inputs.txt");

	string firstline;
	getline(Infile, firstline);

	Infile >> dx >> dt >> dt_Q >> type >> 
		bank_erosion >> dt_output >> meandering >> input_type >> omegac_star >> b >> MC >> max_threads;

	Infile.close();

	std::cout << "XS Spacing (m): " << dx << " | Time Step (s): " << dt << " | Sed Transport EQ: " << type << "\n";
	std::cout << "Discharge time step (s): " << dt_Q << " | Bank erosion: " << bank_erosion << "\n";

	//Start timing
	uint64 time_init = GetTimeMs64();

	//Get number of days from Qfile
	ifstream Qfile;
	Qfile.open(input_path + "Input Q.txt");

	//Count the number of lines and columns
	int n_lines = 0, n_col = 0;
	string line, temp;
	stringstream ss;

	getline(Qfile, line);
	ss.clear();
	ss << line;
	while (ss >> temp) ++n_col;
	++n_lines;

	while (getline(Qfile, line)) ++n_lines;
	Qfile.close();

	//Set number of days and number of time steps
	int n_days = n_lines;
	int n_ts = n_days * dt_Q / dt; //number of time steps for simulation
	int n_nodes = n_col; //number of nodes/reaches

	//Link matrix - specifies how reaches connect, network geometry
	vector<vector<int>> link(2, vector<int>(n_nodes));

	ifstream link_file;
	link_file.open(input_path + "Input link.txt");

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < n_nodes; ++j) {
			link_file >> link[i][j];
		}
	}
	link_file.close();

	//Get reach lengths and number of cross sections
	vector<double> length(n_nodes);
	vector<int> n_xs(n_nodes + 1);
	vector<vector<double>> dx_array(n_nodes, vector<double>(1));
	int max_xs = 0;

	//Initial input function
	input_initial(input_type, input_path, n_nodes, length,
		n_xs, dx_array, max_xs, dx);

	//This set of inputs aren't ever changed
	//Read in Q data from file
	vector<vector<double>> Q(n_days, vector<double>(n_nodes));

	Qfile.open(input_path + "Input Q.txt");
	for (int i = 0; i < n_days; ++i) {
		for (int j = 0; j < n_nodes; ++j) {
			Qfile >> Q[i][j];
		}
	}
	Qfile.close();

	//Grain sizes
	vector<double> Ds;
	ifstream Ds_file;
	Ds_file.open(input_path + "Input Ds.txt");
	
	//count number of lines
	n_lines = 0;
	while (getline(Ds_file, line)) ++n_lines;
	Ds_file.close();

	n_dclass = n_lines;
	Ds.resize(n_dclass);

	Ds_file.open(input_path + "Input Ds.txt");
	double Ds_in;
	int i = 0;
	while (Ds_file >> Ds_in) {
		Ds[i] = Ds_in;
		i = i + 1;
	}
	
	Ds_file.close();

	//Sediment supply
	vector<vector<double>> sed_supply(n_days, vector<double>(n_nodes));

	ifstream sedsupply_file;
	sedsupply_file.open(input_path + "Input sed supply.txt");

	for (int i = 0; i < n_days; ++i) {
		for (int j = 0; j < n_nodes; ++j) {
			sedsupply_file >> sed_supply[i][j];
		}
	}
	sedsupply_file.close();

	//Read in Manning's n values
	vector<double> n_chnl(n_nodes);
	vector<double> n_fp(n_nodes);

	ifstream n_file(input_path + "Input n values.txt");
	for (int i = 0; i < n_nodes; ++i) {
		n_file >> n_chnl[i];
		n_file >> n_fp[i];
	}
	n_file.close();

	//Get all remaining inputs
	vector<vector<double>> bed_z(n_nodes + 1, vector<double>(max_xs));
	vector<vector<double>> bottom_width(n_nodes, vector<double>(max_xs));
	vector<vector<double>> height_RB(n_nodes, vector<double>(max_xs));
	vector<vector<double>> toe_height_RB(n_nodes, vector<double>(max_xs));
	vector<vector<double>> angle_RB(n_nodes, vector<double>(max_xs));
	vector<vector<double>> toe_angle_RB(n_nodes, vector<double>(max_xs));
	vector<vector<double>> height_LB(n_nodes, vector<double>(max_xs));
	vector<vector<double>> toe_height_LB(n_nodes, vector<double>(max_xs));
	vector<vector<double>> angle_LB(n_nodes, vector<double>(max_xs));
	vector<vector<double>> toe_angle_LB(n_nodes, vector<double>(max_xs));
	vector<double> tau_c(n_nodes);
	vector<double> bank_k(n_nodes);
	vector<double> cohesion(n_nodes);
	vector<double> phi(n_nodes);
	vector<double> weight(n_nodes);
	vector<double> cohesion_toe(n_nodes);
	vector<double> phi_toe(n_nodes);
	vector<double> weight_toe(n_nodes);
	vector<double> p_conc(n_nodes);
	vector<double> bed_p_conc(n_nodes);
	vector<double> bank_bedload_prop(n_nodes);
	vector<double> bed_bedload_prop(n_nodes);
	vector<vector<double>> fp_angle(n_nodes, vector<double>(max_xs));
	vector<vector<double>> fp_width_R(n_nodes, vector<double>(max_xs));
	vector<vector<double>> fp_width_L(n_nodes, vector<double>(max_xs));
	vector<vector<double>> bed_tau_c(n_nodes, vector<double>(max_xs));
	vector<vector<double>> cohesive_z(n_nodes, vector<double>(max_xs));
	vector<vector<double>> top_width(n_nodes, vector<double>(max_xs));
	vector<vector<double>> n_bends(n_nodes, vector<double>(max_xs));
	vector<double> sinuosity(n_nodes);
	vector<vector<double>> Rc(n_nodes, vector<double>(max_xs));
	vector<vector<double>> LB_x(n_nodes, vector<double>(max_xs));
	vector<vector<double>> knick_height(n_nodes, vector<double>(max_xs));
	vector<vector<double>> knick_z(n_nodes, vector<double>(max_xs));
	vector<vector<double>> knick_kd(n_nodes, vector<double>(max_xs));
	vector<vector<double>> knick_x(n_nodes, vector<double>(max_xs));
		
	typedef boost::multi_array<double, 3 > array;
	array ps(boost::extents[n_nodes][max_xs][n_dclass]);
	array fs(boost::extents[n_nodes][max_xs][n_dclass]);

	//Get rest of inputs
	input_rest(input_path, input_type, n_nodes, max_xs,
		n_dclass, link, length, n_xs,
		dx_array, bed_z, bottom_width,
		height_RB, toe_height_RB, angle_RB,
		toe_angle_RB, height_LB, toe_height_LB,
		angle_LB, toe_angle_LB, tau_c,
		bank_k, cohesion, phi, weight,
		cohesion_toe, phi_toe, weight_toe, p_conc, bed_p_conc,
		bank_bedload_prop, bed_bedload_prop,
		fp_angle, fp_width_R, fp_width_L, bed_tau_c, cohesive_z,
		top_width, ps, fs, n_bends, sinuosity, Rc, LB_x,
		knick_height, knick_kd, knick_x, knick_z);

	//Factors for adjusting empirical equations in model (sed transport equations, tau_w = f(stream power), 
	//tau_bed = f(stream power), and k = f(tau_w) from Simon et al. 2010)
	double transport_factor = 1, fluvial_factor = 1, cohesive_factor = 1, k_factor = 1;
	//If factor input file, adjust values
	ifstream factor_file(input_path + "Input empirical factors.txt");
	if (factor_file.is_open()) {
		factor_file >> transport_factor >> fluvial_factor >> cohesive_factor >> k_factor;
	}
	factor_file.close();

	if (transport_factor != 1 | fluvial_factor != 1 | cohesive_factor != 1 | k_factor != 1) {
		cout << "You have adjusted the model empirical equations by the following factors:\n";
		cout << "Sediment transport: " << transport_factor << " Fluvial erosion: " << fluvial_factor << "\n";
		cout << "Cohesive bed erosion: " << cohesive_factor << " Soil erodibility: " << k_factor << "\n";
	}

		//If running single model iteration (no Monte Carlo)
		if (MC == 0) {
			//Create output text files
			ofstream zout_file(input_path + "Output z.txt");
			ofstream width_file(input_path + "Output width.txt");
			ofstream dxout_file(input_path + "Output dx.txt");
			ofstream D50_file(input_path + "Output D50.txt");
			ofstream loading_file(input_path + "Output bank loading.txt");
			ofstream geom_file(input_path + "Output XS geometry.txt");

			std::cout << "You are running the normal version of the model.\n";
			//Call main model function
			incision_model(dt, Q, bed_z, link, Ds, bottom_width, top_width, length, sed_supply, lambda, type, nu,
				tau_c, bank_k, dt_output, zout_file, b, width_file, height_RB, height_LB, angle_RB, angle_LB, toe_height_RB,
				toe_height_LB, toe_angle_RB, toe_angle_LB, cohesion, phi, weight, cohesion_toe, phi_toe, weight_toe,
				p_conc, bed_p_conc, bank_erosion, n_xs, dx_array, ps, fs, Rc, sinuosity, 
				n_bends, meandering, dxout_file, cohesive_z, bed_tau_c, omegac_star, D50_file, MC, 
				n_chnl, n_fp, fp_angle, fp_width_R, fp_width_L, bank_bedload_prop, bed_bedload_prop, input_path, loading_file, LB_x,
				geom_file, knick_height, knick_kd, knick_x, knick_z, transport_factor, fluvial_factor, cohesive_factor, 
				k_factor, dt_Q);

			//zout_file.close();
			//width_file.close();
			//D50_file.close();
			//dxout_file.close();
			//geom_file.close();
		}
		//Monte Carlo simulations
		else {
			std::cout << "You are running the Monte Carlo version of the model.\n";

			//For Monte Carlo simulation, determine which inputs are varied
			ifstream MC_file(input_path + "MC Inputs.txt");
			int n_lines = 0;
			while (getline(MC_file, line)) ++n_lines;
			MC_file.close();

			vector<string> MC_inputs(n_lines);
			MC_file.open(input_path + "MC Inputs.txt");
			for (int j = 0; j < n_lines; ++j) {
				MC_file >> MC_inputs[j];
			}
			MC_file.close();

			//Find out which outputs to print
			ifstream MC_outputs_file(input_path + "MC Outputs.txt");
			n_lines = 0;
			while (getline(MC_outputs_file, line)) ++n_lines;
			MC_outputs_file.close();

			vector<string> MC_outputs(n_lines);
			MC_outputs_file.open(input_path + "MC Outputs.txt");
			for (int j = 0; j < n_lines; ++j) {
				MC_outputs_file >> MC_outputs[j];
			}
			MC_outputs_file.close();

			//Get inputs if listed in MC Inputs file
			int n_MC = MC;

			//First need to create all the variables, regardless of whether they will be used
			vector<vector<double>> width_MC(n_nodes, vector<double>(n_MC));
			array ps_MC(boost::extents[n_nodes][n_MC][n_dclass]);
			vector<vector<double>> chnl_n_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> fp_n_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> fp_LB_width_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> fp_RB_width_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> fp_angle_MC(n_nodes, vector<double>(n_MC));
			vector<double> b_MC(n_MC);
			vector<double> omega_star_MC(n_MC);
			vector<vector<double>> tau_c_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> bank_k_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> cohesion_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> phi_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> weight_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> cohesion_toe_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> phi_toe_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> weight_toe_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> bank_p_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> bed_p_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> bank_bedload_prop_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> bed_bedload_prop_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> height_RB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> toe_height_RB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> angle_RB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> toe_angle_RB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> height_LB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> toe_height_LB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> angle_LB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> toe_angle_LB_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> empirical_factors(4, vector<double>(n_MC));
			vector<vector<double>> Ds_MC(n_dclass, vector<double>(n_MC));
			vector<vector<double>> sinuosity_MC(n_nodes, vector<double>(n_MC));
			vector<vector<double>> Rc_MC(n_nodes, vector<double>(n_MC));

			//Variable width
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "width") != MC_inputs.end()) {
				ifstream w_MC_file(input_path + "Input width MC.txt");
				for (int j = 0; j < n_MC; ++j) {
					for (int i = 0; i < n_nodes; ++i) {
						w_MC_file >> width_MC[i][j];
					}
				}
				w_MC_file.close();
			}
			
			//Variable Ds
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "Ds") != MC_inputs.end()) {
				ifstream Ds_MC_file(input_path + "Input Ds MC.txt");
				for (int j = 0; j < n_MC; ++j) {
						for (int k = 0; k < n_dclass; ++k) {
							Ds_MC_file >> Ds_MC[k][j];
						}
				}
				Ds_MC_file.close();
			}

			//Variable GSD
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "ps") != MC_inputs.end()) {
				ifstream ps_MC_file(input_path + "Input ps MC.txt");
				for (int j = 0; j < n_MC; ++j) {
					for (int i = 0; i < n_nodes; ++i) {
						for (int k = 0; k < n_dclass; ++k) {
							ps_MC_file >> ps_MC[i][j][k];
						}
					}
				}
				ps_MC_file.close();
			}
			
			//Variable Manning's n
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "n") != MC_inputs.end()) {
				ifstream n_file(input_path + "Input n MC.txt");
				for (int j = 0; j < n_MC; ++j) {
					for (int i = 0; i < n_nodes; ++i) {
						n_file >> chnl_n_MC[i][j];
						n_file >> fp_n_MC[i][j];
					}
				}
				n_file.close();
			}

			//Variable fp geometry
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "fp") != MC_inputs.end()) {
				ifstream fp_file(input_path + "Input fp MC.txt");
				for (int j = 0; j < n_MC; ++j) {
					for (int i = 0; i < n_nodes; ++i) {
						fp_file >> fp_angle_MC[i][j];
						fp_angle_MC[i][j] *= M_PI / 180;
						fp_file >> fp_LB_width_MC[i][j];
						fp_file >> fp_RB_width_MC[i][j];
					}
				}
				fp_file.close();
			}

			//Variable bank properties
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "bank_prop") != MC_inputs.end()) {
				ifstream bank_prop_file(input_path + "Input bank prop MC.txt");
				for (int j = 0; j < n_MC; ++j) {
					for (int i = 0; i < n_nodes; ++i) {
						bank_prop_file >> tau_c_MC[i][j];
						bank_prop_file >> bank_k_MC[i][j];
						bank_prop_file >> cohesion_MC[i][j];
						bank_prop_file >> cohesion_toe_MC[i][j];
						bank_prop_file >> phi_MC[i][j];
						bank_prop_file >> phi_toe_MC[i][j];
						bank_prop_file >> weight_MC[i][j];
						bank_prop_file >> weight_toe_MC[i][j];
						bank_prop_file >> bank_p_MC[i][j];
						bank_prop_file >> bed_p_MC[i][j];
						bank_prop_file >> bank_bedload_prop_MC[i][j];
						bank_prop_file >> bed_bedload_prop_MC[i][j];
					}
				}
				bank_prop_file.close();
			}

			//Variable bank geometry
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "geom") != MC_inputs.end()) {
				ifstream RB_geom_file(input_path + "Input RB geometry MC.txt");
				ifstream LB_geom_file(input_path + "Input LB geometry MC.txt");
				for (int j = 0; j < n_MC; ++j) {
					for (int i = 0; i < n_nodes; ++i) {
						RB_geom_file >> height_RB_MC[i][j];
						RB_geom_file >> toe_height_RB_MC[i][j];
						RB_geom_file >> angle_RB_MC[i][j];
						RB_geom_file >> toe_angle_RB_MC[i][j];

						LB_geom_file >> height_LB_MC[i][j];
						LB_geom_file >> toe_height_LB_MC[i][j];
						LB_geom_file >> angle_LB_MC[i][j];
						LB_geom_file >> toe_angle_LB_MC[i][j];
					}
				}
				RB_geom_file.close();
				LB_geom_file.close();
			}

			//Get omega_c_star and b for hiding function
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "hiding") != MC_inputs.end()) {
				ifstream hiding_file(input_path + "Input hiding MC.txt");

				for (int i = 0; i < n_MC; ++i) {
					hiding_file >> omega_star_MC[i];
					hiding_file >> b_MC[i];
				}
				hiding_file.close();
			}

			//Variable meandering
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "meandering") != MC_inputs.end()) {
				ifstream meandering_file(input_path + "Input meandering MC.txt");
				for (int j = 0; j < n_MC; ++j) {
					for (int i = 0; i < n_nodes; ++i) {
						meandering_file >> sinuosity_MC[i][j];
						meandering_file >> Rc_MC[i][j];
					}
				}
				meandering_file.close();
			}

			//Variable empirical factors
			if (std::find(MC_inputs.begin(), MC_inputs.end(), "factors") != MC_inputs.end()) {
				ifstream factor_file(input_path + "Input empirical factors MC.txt");

				for (int i = 0; i < n_MC; ++i) {
					factor_file >> empirical_factors[0][i] >> empirical_factors[1][i] >> empirical_factors[2][i] >>
						empirical_factors[3][i];
				}
				factor_file.close();
			}

			int iteration = 0;
			//Set up progress bar
			std::cout << "\n\n" << "0%____________________100%\n" << "  ";
//Parallelize model iterations to increase speed. Need to set all variables to private so they aren't
//changed simultaneously by the different processesors
			if (max_threads != 0) {
				omp_set_num_threads(max_threads);
			}
			//int n_processors = omp_get_max_threads();
			//if (n_processors > 35) {
			//	omp_set_num_threads(35);
			//}
			//dt_output, bank_erosion, meandering, MC, transport_factor, fluvial_factor, cohesive_factor, k_factor

#pragma omp parallel for default(none) \
	firstprivate(bed_z, bottom_width, top_width, length, tau_c, bank_k, b, height_RB, \
	height_LB, angle_RB, angle_LB, toe_height_RB, toe_height_LB, toe_angle_RB, toe_angle_LB, cohesion, \
	phi, weight, cohesion_toe, phi_toe, weight_toe, dx_array, ps, fs, Rc, sinuosity, n_bends, cohesive_z, \
	bed_tau_c, omegac_star, n_chnl,	n_fp, fp_angle, fp_width_R, fp_width_L, b_MC, omega_star_MC, width_MC, fp_angle_MC, \
	fp_LB_width_MC, fp_RB_width_MC, chnl_n_MC, fp_n_MC, ps_MC, dt, Q, link, Ds, sed_supply, lambda, type, nu, n_xs, p_conc, LB_x, \
	bed_p_conc, bank_bedload_prop, bed_bedload_prop, knick_height, knick_kd, knick_x, knick_z, MC_outputs, MC_inputs, \
	input_path, Ds_MC, sinuosity_MC, Rc_MC, dt_Q)
				for (int i_MC = 0; i_MC < n_MC; ++i_MC) {
					//std::cout << "Thread " << omp_get_thread_num() << " Iteration " << i + 1 << "\n";
					iteration += 1;
					//Reset all inputs
					//Reset n_xs
					for (int i = 0; i <= n_nodes; ++i) {
						n_xs[i] = 0;
					}

					for (int j = 0; j < n_nodes; ++j) {
						for (int k = 0; k < max_xs; ++k) {
							knick_height[j][k] = 0;
							knick_kd[j][k] = 0;
							knick_x[j][k] = 0;
							knick_z[j][k] = 0;
						}
					}

					input_initial(input_type, input_path, n_nodes, length, n_xs, dx_array, max_xs, dx);
					input_rest(input_path, input_type, n_nodes, max_xs,
						n_dclass, link, length, n_xs,
						dx_array, bed_z, bottom_width,
						height_RB, toe_height_RB, angle_RB,
						toe_angle_RB, height_LB, toe_height_LB,
						angle_LB, toe_angle_LB, tau_c,
						bank_k, cohesion, phi, weight,
						cohesion_toe, phi_toe, weight_toe, p_conc, bed_p_conc,
						bank_bedload_prop, bed_bedload_prop,
						fp_angle, fp_width_R, fp_width_L, bed_tau_c, cohesive_z,
						top_width, ps, fs, n_bends, sinuosity, Rc, LB_x,
						knick_height, knick_kd, knick_x, knick_z);
					
					//Hiding function variables
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "hiding") != MC_inputs.end()) {
						b = b_MC[i_MC];
						omegac_star = omega_star_MC[i_MC];
					}

					//Update bottom_width and ps, also Manning's n and fp geometry
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "width") != MC_inputs.end()) {
						for (int j = 0; j < n_nodes; ++j) {
							bottom_width[j][0] = width_MC[j][i_MC];

							//populate rest of xs
							for (int k = 0; k < n_xs[j]; ++k) {
								bottom_width[j][k] = bottom_width[j][0];
							}
						}
					}

					//fp geometry
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "fp") != MC_inputs.end()) {
						for (int j = 0; j < n_nodes; ++j) {
							fp_angle[j][0] = fp_angle_MC[j][i_MC];
							fp_width_R[j][0] = fp_RB_width_MC[j][i_MC];
							fp_width_L[j][0] = fp_LB_width_MC[j][i_MC];

							for (int k = 1; k < n_xs[j]; ++k) {
								fp_angle[j][k] = fp_angle[j][0];
								fp_width_R[j][k] = fp_width_R[j][0];
								fp_width_L[j][k] = fp_width_L[j][0];
							}
						}
					}

					//Manning's n
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "n") != MC_inputs.end()) {
						for (int j = 0; j < n_nodes; ++j) {
							n_chnl[j] = chnl_n_MC[j][i_MC];
							n_fp[j] = fp_n_MC[j][i_MC];
						}
					}

					//Update bank properties
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "bank_prop") != MC_inputs.end()) {
						for (int j = 0; j < n_nodes; ++j) {
							tau_c[j] = tau_c_MC[j][i_MC];
							bank_k[j] = bank_k_MC[j][i_MC];
							cohesion[j] = cohesion_MC[j][i_MC];
							phi[j] = phi_MC[j][i_MC] * M_PI / 180.0;
							weight[j] = weight_MC[j][i_MC];
							cohesion_toe[j] = cohesion_toe_MC[j][i_MC];
							phi_toe[j] = phi_toe_MC[j][i_MC] * M_PI / 180.0;
							weight_toe[j] = weight_toe_MC[j][i_MC];
							p_conc[j] = bank_p_MC[j][i_MC];
							bed_p_conc[j] = bed_p_MC[j][i_MC];
							bank_bedload_prop[j] = bank_bedload_prop_MC[j][i_MC];
							bed_bedload_prop[j] = bed_bedload_prop_MC[j][i_MC];
						}
					}

					//Grain size distribution
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "ps") != MC_inputs.end()) {
						for (int j = 0; j < n_nodes; j++) {
							for (int m = 0; m < n_dclass; m++) {
								ps[j][0][m] = ps_MC[j][i_MC][m];
								//Populate all XS in reach
								for (int k = 0; k < n_xs[j]; ++k) {
									ps[j][k][m] = ps[j][0][m];
									fs[j][k][m] = ps[j][0][m];
								}
							}
						}
					}

					//Grain sizes
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "Ds") != MC_inputs.end()) {
						for (int m = 0; m < n_dclass; m++) {
							Ds[m] = Ds_MC[m][i_MC];
						}
					}

					//Update bank geometry
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "geom") != MC_inputs.end()) {
						for (int j = 0; j < n_nodes; ++j) {
							height_RB[j][0] = height_RB_MC[j][i_MC];
							toe_height_RB[j][0] = toe_height_RB_MC[j][i_MC];
							//Convert angles to radians
							angle_RB[j][0] = angle_RB_MC[j][i_MC] * M_PI / 180;
							toe_angle_RB[j][0] = toe_angle_RB_MC[j][i_MC] * M_PI / 180;

							height_LB[j][0] = height_LB_MC[j][i_MC];
							toe_height_LB[j][0] = toe_height_LB_MC[j][i_MC];
							//Convert angles to radians
							angle_LB[j][0] = angle_LB_MC[j][i_MC] * M_PI / 180;
							toe_angle_LB[j][0] = toe_angle_LB_MC[j][i_MC] * M_PI / 180;

							//populate rest of xs
							for (int k = 0; k < n_xs[j]; ++k) {
								height_RB[j][k] = height_RB[j][0];
								toe_height_RB[j][k] = toe_height_RB[j][0];
								angle_RB[j][k] = angle_RB[j][0];
								toe_angle_RB[j][k] = toe_angle_RB[j][0];

								height_LB[j][k] = height_LB[j][0];
								toe_height_LB[j][k] = toe_height_LB[j][0];
								angle_LB[j][k] = angle_LB[j][0];
								toe_angle_LB[j][k] = toe_angle_LB[j][0];
							}
						}
					}

					//Calculate top width based on bottom width, bank height, and bank angle
					vector<vector<double>> top_width(n_nodes, vector<double>(max_xs));
					for (int i = 0; i < n_nodes; ++i) {
						for (int j = 0; j < n_xs[i]; ++j) {
							top_width[i][j] = bottom_width[i][j] + toe_height_RB[i][j] / tan(toe_angle_RB[i][j])
								+ (height_RB[i][j] - toe_height_RB[i][j]) / tan(angle_RB[i][j])
								+ toe_height_LB[i][j] / tan(toe_angle_LB[i][j])
								+ (height_LB[i][j] - toe_height_LB[i][j]) / tan(angle_LB[i][j]);
						}
					}

					//Variable meandering
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "meandering") != MC_inputs.end()) {
						for (int j = 0; j < n_nodes; ++j) {
							sinuosity[j] = sinuosity_MC[j][i_MC];
							Rc[j][0] = Rc_MC[j][i_MC];

							//Calculate num_bends and fill in Rc values for rest of XS
							typedef pair<double, double> Result;
							for (int k = 0; k < n_xs[j]; ++k) {
								Result n_bends_calc = boost::math::tools::brent_find_minima(
									bind(calc_n_bends, placeholders::_1, dx_array[j][k], Rc[j][0], sinuosity[j]), 0.0, 
									1000.0, 20);

								n_bends[j][k] = n_bends_calc.first;
								Rc[j][k] = Rc[j][0];
							}
						}
	
					}

					//Variable empirical factors
					if (std::find(MC_inputs.begin(), MC_inputs.end(), "factors") != MC_inputs.end()) {
						transport_factor = empirical_factors[0][i_MC]; 
						fluvial_factor = empirical_factors[1][i_MC];
						cohesive_factor = empirical_factors[2][i_MC];
						k_factor = empirical_factors[3][i_MC];
					}

					//Create output files
					ofstream zout_file;
					ofstream D50_file;
					ofstream loading_file;
					ofstream width_file;

					int open_fail = 0;
					if (std::find(MC_outputs.begin(), MC_outputs.end(), "z") != MC_outputs.end()) {
						zout_file.open(input_path + "MC Outputs/Output z" + std::to_string(i_MC) + ".txt");
						//ifstream test_file(input_path + "MC Outputs/Output z" + std::to_string(i_MC) + ".txt");
						if (!zout_file.is_open()) {
							open_fail += 1;
							//cout << "File " << i_MC << " didn't open.\n";
							//cout << input_path + "MC Outputs/Output z" + std::to_string(i_MC) + ".txt\n";
							//cout << zout_file << "\n";
							//zout_file.open(input_path + "MC Outputs/Output z" + std::to_string(i_MC) + ".txt");
						}
					}
					else {
						zout_file.open(input_path + "Output z.txt");
					}

					if (std::find(MC_outputs.begin(), MC_outputs.end(), "D50") != MC_outputs.end()) {
						D50_file.open(input_path + "MC Outputs/Output D50" + std::to_string(i_MC) + ".txt");
						if (!D50_file.is_open()) {
							open_fail += 1;
						}
					}
					else {
						D50_file.open(input_path + "Output D50.txt");
					}

					if (std::find(MC_outputs.begin(), MC_outputs.end(), "loading") != MC_outputs.end()) {
						loading_file.open(input_path + "MC Outputs/Output bank loading" + std::to_string(i_MC) + ".txt");
						if (!loading_file.is_open()) {
							open_fail += 1;
						}
					}
					else {
						loading_file.open(input_path + "Output bank loading.txt");
					}

					if (std::find(MC_outputs.begin(), MC_outputs.end(), "width") != MC_outputs.end()) {
						width_file.open(input_path + "MC Outputs/Output width" + std::to_string(i_MC) + ".txt");
						if (!width_file.is_open()) {
							open_fail += 1;
						}
					}
					else {
						width_file.open(input_path + "Output width.txt");
					}
					
					if (open_fail > 0) {
						std::cout << "\nOne or more output files failed to open. Try reducing 'max_threads'.\n";
					}

					//I don't care about these outputs right now
					ofstream dxout_file(input_path + "Output dx.txt");
					ofstream geom_file(input_path + "Output XS geometry.txt");

					//Call model
					incision_model(dt, Q, bed_z, link, Ds, bottom_width, top_width, length, sed_supply, lambda, type, nu,
						tau_c, bank_k, dt_output, zout_file, b, width_file, height_RB, height_LB, angle_RB, angle_LB, toe_height_RB,
						toe_height_LB, toe_angle_RB, toe_angle_LB, cohesion, phi, weight, cohesion_toe, phi_toe, weight_toe,
						p_conc, bed_p_conc, bank_erosion, n_xs, dx_array, ps, fs, Rc, sinuosity,
						n_bends, meandering, dxout_file, cohesive_z, bed_tau_c, omegac_star, D50_file, MC,
						n_chnl, n_fp, fp_angle, fp_width_R, fp_width_L, bank_bedload_prop, bed_bedload_prop, input_path, loading_file, LB_x,
						geom_file, knick_height, knick_kd, knick_x, knick_z, transport_factor, fluvial_factor, 
						cohesive_factor, k_factor, dt_Q);

					//Print progress
					if ((iteration / max(1, (n_MC / 20))) != ((iteration - 1) / max(1, (n_MC / 20)))) {
						std::cout << "=";
					}
				}
			
		}

	//Print final summary file
	print_summary_file(input_path, time_init, dx, dt, type, bank_erosion, MC);
}

void print_summary_file(string& input_path, uint64& time_init, int& dx, int& dt, string& type, string& bank_erosion,
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
	out_file << "WARNING MESSAGES\n";
	//cout << "\n\n" << "Maximum bed elevation change in time step (abs): " << max_dz << "\n\n";

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
	//ShellExecuteW(NULL, "open", input_path + "Model Output Summary.txt", "", input_path, SW_SHOW);
}