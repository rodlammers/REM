#include "stdafx.h"

using namespace std;

void input_initial(string& input_type, const string& input_path, int& n_nodes,
	vector<double>& length, vector<int>& n_xs, 
	vector<vector<double>>& dx_array, int& max_xs, int& dx) {
	
	ifstream length_file(input_path + "Input length.txt");

	////Reset n_xs
	//for (int i = 0; i <= n_nodes; ++i) {
	//	n_xs[i] = 0;
	//}
	
	if (input_type == "reach") {
		//if (meandering == 1) {
		//	ifstream sinuosity_file(input_path + "Input sinuosity.txt");
		//	ifstream Rc_file(input_path + "Input Rc.txt");

		//	for (int i = 0; i < n_nodes; ++i) {
		//		sinuosity_file >> sinuosity[i];
		//		Rc_file >> Rc[i][0];
		//		length_file >> length[i];

		//		//Calculate num_bends and num_xs
		//		typedef pair<double, double> Result;
		//		Result n_bends_calc = boost::math::tools::brent_find_minima(
		//			bind(calc_n_bends, placeholders::_1, length[i], Rc[i][0], sinuosity[i]), 1.0, 1000.0, 20);

		//		n_bends[i] = n_bends_calc.first;
		//		n_xs[i] = max(floor(n_bends[i]) * 2.0, 1.0); //Number of XS is twice the number of bends plus one
		//		if (fmod(n_bends[i], floor(n_bends[i])) < 5) {
		//			n_xs[i] += 1;
		//		}
		//		else
		//		{
		//			n_xs[i] += 2;
		//		}
		//		dx_array[i][0] = length[i] / (n_bends[i] * 2);
		//		if (n_xs[i] > max_xs) { max_xs = n_xs[i]; }

		//	}

		//	sinuosity_file.close();
		//	Rc_file.close();

		//	ofstream Rc_out_file(input_path + "Output Rc.txt");
		//	//Resize Rc vector and add values for remaining XS's
		//	for (int i = 0; i < n_nodes; ++i) {
		//		Rc[i].resize(max_xs);
		//		dx_array[i].resize(max_xs);
		//		Rc_out_file << Rc[i][0] << "	";
		//		for (int j = 1; j < n_xs[i]; ++j) {
		//			if (j % 2 == 0) {
		//				Rc[i][j] = 0;
		//			}
		//			else {
		//				Rc[i][j] = Rc[i][0];
		//			}
		//			dx_array[i][j] = dx_array[i][0];
		//			Rc_out_file << Rc[i][j] << "	";
		//		}
		//		//Reset initial Rc to zero
		//		Rc[i][0] = 0;
		//		Rc_out_file << "\n";
		//	}
		//	Rc_out_file.close();
		//}
		//else {
			for (int i = 0; i < n_nodes; ++i) {
				length_file >> length[i];

				n_xs[i] = max(round(length[i] / dx) + 1, 2.0);
				//n_xs[i] = max(round(length[i] / dx), 1);

				dx_array[i][0] = length[i] / (n_xs[i] - 1);
				n_xs[i] -= 1;
				if (n_xs[i] > max_xs) { max_xs = n_xs[i]; }
			}

			//Populate remainder of dx_array
			for (int i = 0; i < n_nodes; ++i) {
				dx_array[i].resize(max_xs);
				for (int j = 1; j < n_xs[i]; ++j) {
					dx_array[i][j] = dx_array[i][0];
				}
			}

		//}

		length_file.close();
	}
	else {
		ifstream z_file(input_path + "Input z.txt");

		int n_lines = 0;
		string line;
		while (getline(z_file, line)) ++n_lines;
		max_xs = n_lines;

		z_file.close();
	}

	//Set last n_xs to one
	n_xs[n_nodes] = 1;
	
}


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
	vector<vector<double>>& knick_x, vector<vector<double>>& knick_z) {

	////Reset n_xs
	//for (int i = 0; i <= n_nodes; ++i) {
	//	n_xs[i] = 0;
	//}

	ifstream z_file(input_path + "Input z.txt");
	if (input_type == "reach") {
		double z;
		int i = 0;
		while (z_file >> z) {
			bed_z[i][0] = z; //Gets elevation of each node
			i = i + 1;
		}

		//Now populate rest of XS elevations
		for (int i = 0; i < n_nodes; ++i) {
			//Find DS node index
			int DS_index = find(link, i + 1);
			double slope = (bed_z[i][0] - bed_z[DS_index][0]) / length[i];
			for (int j = 1; j < n_xs[i]; ++j) {
				bed_z[i][j] = bed_z[i][j - 1] - slope * dx_array[i][j];
			}
		}

		//Read in cohesive bed information
		ifstream cohesive_bed_file(input_path + "Input bed cohesive.txt");
		for (int i = 0; i < n_nodes; ++i) {
			double cohesive_depth;
			cohesive_bed_file >> bed_tau_c[i][0];
			cohesive_bed_file >> cohesive_depth;

			//Convert cohesive depth to elevation, for all XS in a reach
			for (int j = 0; j < n_xs[i]; ++j) {
				cohesive_z[i][j] = bed_z[i][j] - cohesive_depth;
				bed_tau_c[i][j] = bed_tau_c[i][0];
			}
		}
		cohesive_bed_file.close();
	}
	else {
		for (int j = 0; j < max_xs; ++j) {
			for (int i = 0; i < n_nodes + 1; ++i) {
				z_file >> bed_z[i][j];

				if (bed_z[i][j] != 0) n_xs[i] += 1;
			}
		}

		//Read in dx values from that file
		ifstream dx_file(input_path + "Input dx.txt");

		for (int j = 0; j < max_xs; ++j) {
			for (int i = 0; i < n_nodes; ++i) {
				if (j == 0) {
					dx_array[i].resize(max_xs);
				}
				dx_file >> dx_array[i][j];
			}
		}

		dx_file.close();

		//Read in cohesive bed information
		ifstream cohesive_bed_file(input_path + "Input bed cohesive.txt");
		for (int j = 0; j < max_xs; ++j) {
			double cohesive_depth;
			//Convert cohesive depth to elevation, for all XS in a reach
			for (int i = 0; i < n_nodes; ++i) {
				cohesive_bed_file >> bed_tau_c[i][j];
				cohesive_bed_file >> cohesive_depth;
				cohesive_z[i][j] = bed_z[i][j] - cohesive_depth;
			}
		}
		cohesive_bed_file.close();
	}
	z_file.close();

	//Set last n_xs to one
	n_xs[n_nodes] = 1;

	//Input sinuosity
	ifstream meandering_file(input_path + "Input meandering.txt");

	for (int i = 0; i < n_nodes; ++i) {
		meandering_file >> sinuosity[i];
		meandering_file >> Rc[i][0];

		//Calculate num_bends and fill in Rc values for rest of XS
		typedef pair<double, double> Result;
		for (int j = 0; j < n_xs[i]; ++j) {
			Result n_bends_calc = boost::math::tools::brent_find_minima(
				bind(calc_n_bends, placeholders::_1, dx_array[i][j], Rc[i][0], sinuosity[i]), 0.0, 1000.0, 20);

			n_bends[i][j] = n_bends_calc.first;

			Rc[i][j] = Rc[i][0];
		}
	}

	meandering_file.close();

	//Channel width input (m)
	//Width from file
	ifstream w_file;
	w_file.open(input_path + "Input width.txt");

	for (int i = 0; i < n_nodes; ++i) {
		w_file >> bottom_width[i][0];
		//Populate rest of XS
		for (int j = 1; j < n_xs[i]; ++j) {
			bottom_width[i][j] = bottom_width[i][0];
		}
	}
	w_file.close();

	//Bank geometry - RB
	ifstream RB_file(input_path + "Input RB geometry.txt");

	for (int i = 0; i < n_nodes; ++i) {
		RB_file >> height_RB[i][0];
		RB_file >> toe_height_RB[i][0];
		RB_file >> angle_RB[i][0];
		RB_file >> toe_angle_RB[i][0];

		//Convert to radians
		angle_RB[i][0] = angle_RB[i][0] * M_PI / 180;
		toe_angle_RB[i][0] = toe_angle_RB[i][0] * M_PI / 180;

		//Populate rest of XS
		for (int j = 1; j < n_xs[i]; ++j) {
			height_RB[i][j] = height_RB[i][0];
			toe_height_RB[i][j] = toe_height_RB[i][0];
			angle_RB[i][j] = angle_RB[i][0];
			toe_angle_RB[i][j] = toe_angle_RB[i][0];

		}
	}
	RB_file.close();

	//Bank geometry - LB
	ifstream LB_file(input_path + "Input LB geometry.txt");

	for (int i = 0; i < n_nodes; ++i) {
		LB_file >> height_LB[i][0];
		LB_file >> toe_height_LB[i][0];
		LB_file >> angle_LB[i][0];
		LB_file >> toe_angle_LB[i][0];

		//Convert to radians
		angle_LB[i][0] = angle_LB[i][0] * M_PI / 180;
		toe_angle_LB[i][0] = toe_angle_LB[i][0] * M_PI / 180;

		//Populate rest of XS
		for (int j = 1; j < n_xs[i]; ++j) {
			height_LB[i][j] = height_LB[i][0];
			toe_height_LB[i][j] = toe_height_LB[i][0];
			angle_LB[i][j] = angle_LB[i][0];
			toe_angle_LB[i][j] = toe_angle_LB[i][0];

		}
	}
	LB_file.close();

	//LB station - location of the bottom of the left toe point. Tracked throughout the model to
	//allow for XS plotting of outputs.
	ifstream LBx_file(input_path + "Input LB x.txt");
	for (int i = 0; i < n_nodes; ++i) {
		LBx_file >> LB_x[i][0];
		for (int j = 1; j < n_xs[i]; ++j) {
			LB_x[i][j] = LB_x[i][0];
		}
	}
	LBx_file.close();

	//Read in bank soil properties
	ifstream bank_prop_file(input_path + "Input bank prop.txt");
	for (int i = 0; i < n_nodes; ++i) {
		bank_prop_file >> tau_c[i];
		bank_prop_file >> bank_k[i];
		bank_prop_file >> cohesion[i];
		bank_prop_file >> cohesion_toe[i];
		bank_prop_file >> phi[i];
		phi[i] = phi[i] * M_PI / 180.0; //convert to radians
		bank_prop_file >> phi_toe[i];
		phi_toe[i] = phi_toe[i] * M_PI / 180.0; //convert to radians
		bank_prop_file >> weight[i];
		bank_prop_file >> weight_toe[i];
		bank_prop_file >> p_conc[i];
		bank_prop_file >> bed_p_conc[i];
		bank_prop_file >> bank_bedload_prop[i];
		bank_prop_file >> bed_bedload_prop[i];
	}
	bank_prop_file.close();

	//Read in floodplain geometry
	ifstream fp_geom_file(input_path + "Input fp geometry.txt");
	for (int i = 0; i < n_nodes; ++i) {
		fp_geom_file >> fp_angle[i][0];
		fp_angle[i][0] = fp_angle[i][0] * M_PI / 180.0; //convert to radians
		fp_geom_file >> fp_width_L[i][0];
		fp_geom_file >> fp_width_R[i][0];

		for (int j = 0; j < n_xs[i]; ++j) {
			fp_angle[i][j] = fp_angle[i][0];
			fp_width_R[i][j] = fp_width_R[i][0];
			fp_width_L[i][j] = fp_width_L[i][0];
		}
	}
	fp_geom_file.close();

	////Initialize ps arrays
	ifstream ps_file;
	ps_file.open(input_path + "Input ps.txt");

	for (int j = 0; j < n_nodes; j++) {
		for (int m = 0; m < n_dclass; m++) {
			ps_file >> ps[j][0][m];
			//Populate all XS in reach
			for (int k = 0; k < n_xs[j]; ++k) {
				ps[j][k][m] = ps[j][0][m];
				fs[j][k][m] = ps[j][0][m];
			}
		}
	}
	ps_file.close();

	//Calculate top width based on bottom width, bank height, and bank angle
	for (int i = 0; i < n_nodes; ++i) {
		for (int j = 0; j < n_xs[i]; ++j) {
			top_width[i][j] = bottom_width[i][j] + toe_height_RB[i][j] / tan(toe_angle_RB[i][j])
				+ (height_RB[i][j] - toe_height_RB[i][j]) / tan(angle_RB[i][j])
				+ toe_height_LB[i][j] / tan(toe_angle_LB[i][j])
				+ (height_LB[i][j] - toe_height_LB[i][j]) / tan(angle_LB[i][j]);
		}
	}

	//Get knickpoint info
	vector < vector<double>> knick_dist(n_nodes, vector<double>(max_xs));
	ifstream knick_file(input_path + "Input knickpoint.txt");
	int j = 0;
	while (knick_file >> j) {
		//j is reach with knickpoint. Find which XS knickpoint is ds of
		j -= 1;
		double dist;
		knick_file >> dist; //DS distance
		int xs = 0, cum_dist = 0;
		for (int k = 0; k < n_xs[j]; ++k) {
			cum_dist += dx_array[j][k];
			if (cum_dist > dist) {
				xs = k;
				cum_dist -= dx_array[j][k];
				break;
			}
		}
		knick_file >> knick_z[j][xs];
		knick_file >> knick_height[j][xs];
		knick_file >> knick_kd[j][xs];
		knick_x[j][xs] = dist - cum_dist;
		knick_dist[j][xs] = dist;
	}
	knick_file.close();

	if (input_type == "reach") {
		//Now populate rest of XS elevations
		for (int j = 0; j < n_nodes; ++j) {
			//Find DS node index
			int DS_index = find(link, j + 1);
			double slope = (bed_z[j][0] - bed_z[DS_index][0]) / length[j];
			int US_k = 0;
			for (int k = 1; k < n_xs[j]; ++k) {
				double US_elev = bed_z[j][0];
				if (knick_height[j][k - 1] != 0) {
					slope = (US_elev - knick_z[j][k-1]) / knick_dist[j][k-1];
					for (int m = US_k + 1; m < k; ++m) {
						double diff = bed_z[j][m] - (bed_z[j][m - 1] - slope * dx_array[j][m]);
						bed_z[j][m] += diff;
						cohesive_z[j][m] += diff;
					}
					US_elev = knick_z[j][k - 1];
					US_k = k - 1;
				}
			}
			if (US_k == 0 & knick_z[j][n_xs[j] - 1] != 0) {
				US_k = n_xs[j] - 1;
				//slope = (knick_z[j][US_k] - knick_height[j][US_k] - bed_z[DS_index][0]) / (length[j] - knick_dist[j][US_k]);
				slope = (bed_z[j][0] - knick_z[j][US_k]) / (knick_dist[j][US_k]);
				double diff = bed_z[j][n_xs[j] - 1] - (bed_z[DS_index][0] + slope * dx_array[j][n_xs[j] - 1]);
				bed_z[j][n_xs[j] - 1] += diff;
				cohesive_z[j][n_xs[j] - 1] += diff;
				for (int k = n_xs[j] - 2; k > US_k; --k) {
					diff = bed_z[j][k] - (bed_z[j][k + 1] + slope * dx_array[j][k]);
					bed_z[j][k] -= diff;
					cohesive_z[j][k] -= diff;
				}
			}
		}
	}
}