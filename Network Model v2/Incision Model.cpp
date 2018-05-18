#include "stdafx.h"

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
	double& k_factor, int& dt_Q) {

	int n_days = Q.size(); //number of days to simulation
	long n_ts = n_days * dt_Q / dt; //number of time steps for simulation
	int n_nodes = bed_z.size() - 1; //number of nodes (reaches)
	int Q_iter, Q_iter_old = -1; //simulation iteration by discharge (daily, hourly, 15-min, etc)
	const double gamma = 9810; //specific weight of water (N/m^3)
	const double a = -1; //0.5 = Quick scheme; 1 = central difference; -1 = second order upwind
	//int dt_bank = 86400; //seconds in a day
	//int dt_bank = dt;
	int n_xs_max = bed_z[0].size(); //maximum number of XS in a reach
	double max_dz = 0; //maximum bed elevation change in a single time step
	int n_dclass = Ds.size(); //number of grain size classes
	string limit_fun = "Superbee"; //selected limiter function

	vector<double> Ds_log(n_dclass); //log2 of grain sizes
	vector<vector<double>> ds_old(n_nodes, vector<double>(n_xs_max)); //Old active layer thickness (m)
	vector<vector<double>> q(n_nodes, vector<double>(n_xs_max)); //unit discharge (m^3/s/m)
	vector<vector<double>> slope(n_nodes, vector<double>(n_xs_max)); //channel slope (m/m)
	vector<vector<double>> bank_sed(n_nodes, vector<double>(n_xs_max)); //sediment loading rate from bank erosion (m^3/s)
	vector<vector<double>> cohesive_bedload(n_nodes, vector<double>(n_xs_max)); //bedload loading rate from cohesive bed erosion (m^3/s)
	vector<vector<double>> knick_bedload(n_nodes, vector<double>(n_xs_max)); //bedload loading rate from knickpoint erosion (m^3/s)
	vector<vector<double>> bank_sed_mass(n_nodes, vector<double>(n_xs_max)); //mass sediment loading from bank erosion (kg)
	vector<vector<double>> bank_p_mass(n_nodes, vector<double>(n_xs_max)); //mass phosphorus loading from bank erosion (kg)
	vector<double> sed_supply_day(n_nodes); //Vector of only sed supply for that day
	vector<double> sed_supply_output(n_nodes); //Vector of calculated sediment inputs for each reach, to be outputted for use (m^3/m/s)
	vector<vector<double>> bank_tank_RB(n_nodes, vector<double>(n_xs_max)); //"tank" of failed bank material stored at the
																		//bank toe. The bank isn't eroded further until this material is removed.
	vector<vector<double>> bank_tank_LB(n_nodes, vector<double>(n_xs_max));
	vector<vector<double>> cohesive_sed_mass(n_nodes, vector<double>(n_xs_max)); //mass sed loading from cohesive bed erosion (kg)
	vector<vector<double>> knick_sed_mass(n_nodes, vector<double>(n_xs_max)); //mass sed loading from knickpoint migration (kg)
	vector<vector<double>>bed_p_mass(n_nodes, vector<double>(n_xs_max)); //mass phosphorus loading from cohesive bed erosion (kg)

	vector<vector<double>> meander_erosion(n_nodes, vector<double>(n_xs_max));
	vector<vector<double>> avg_width(n_nodes, vector<double>(n_xs_max)); // Average bottom width of XS and immediately DS XS (m)

	double bed_sed_out = 0;  //Mass of bed material load sediment exported from watershed (kg)
	double bed_sed_out_vol = 0; //Volume of bed material load sediment exported from watershed (m^3)
	double bank_sed_vol = 0; //Volume of bank-derived sediment added to bed material load (m^3)
	double bed_sed_in_vol = 0; //Volume of any sediment input (m^3)
	double bank_wash_vol = 0; //Volume of washload from bank erosion (m^3)
	double cohesive_bed_vol = 0; //Volume of washload from cohesive bed erosion (m^3)
	double knick_bed_vol = 0; //Volume of washload from knickpoint migration (m^3)
	double knick_vol_correction = 0; //Volume correction to account for eroded area in XS downstream from where knickpoint
									//originates. This isn't added to washload but is used to balance eroded/deposited volumes.
	double dQ_dx_sum = 0; //error in calculated change in dz (m)
	vector<vector<int>> skip_bank_erosion_L(n_nodes, vector<int>(n_xs_max)); //Whether or not to skip bank erosion for
																			 //that location (if aggraded out to edge of fp.)
																			 //1 = no failure and 2 = no failure or fluvial
	vector<vector<int>> skip_bank_erosion_R(n_nodes, vector<int>(n_xs_max));

	//Output files
	//ofstream qb_file(input_path + "Output Qb.txt");
	//ofstream bank_sed_out_file(input_path + "Output bank_sed.txt");
	ofstream sinuosity_out_file(input_path + "Output sinuosity.txt");
	ofstream geom_file2(input_path + "Output XS geometry all.txt");
	ofstream stream_power_file(input_path + "Output stream power.txt");
	ofstream slope_file(input_path + "Output slope.txt");
	ofstream knick_out_file(input_path + "Output knick locations.txt");
	//ofstream sed_supply_output_file(input_path + "Output sed supply.txt");
	ofstream Rc_file(input_path + "Output Rc.txt");

	ofstream bed_sed_mass_file;
	ofstream sediment_volume_file;

	//Remove error file, if present
	remove((input_path + "MODEL ERRORS.txt").c_str());

	//Print headers for bank loading file
	for (int j = 0; j < n_nodes; ++j) {
		loading_file << "Sed" + to_string(j + 1) << "	";
	}
	for (int j = 0; j < n_nodes; ++j) {
		loading_file << "P" + to_string(j + 1) << "	";
	}
	loading_file << "\n";

	if (MC == 0) {
		bed_sed_mass_file.open(input_path + "Output bed mass.txt");
		sediment_volume_file.open(input_path + "Output sediment vol.txt");

		sediment_volume_file << "Bed_vol_out	Bed_vol_in	Bank_tank_vol	Bank_washload	" <<
			"Cohesive_bed_washload	Knickpoint_washload	Bank_bedload_vol	Knick_correction	Calc_error\n";
	}


	if (MC == 0) {
		//Set up progress bar
		cout << "\n\n" << "0%____________________100%\n" << "  ";
	}

	//Get log of grain sizes
	for (int m = 0; m < n_dclass; ++m) {
		Ds_log[m] = log2(Ds[m]);
	}

	//Add initial eroded knickpoint volume based on remainder of XS
	//Also, if bed tau_c = 999, this is a hardpoint so we turn off bank erosion
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			if (knick_height[j][k] != 0) {
				knick_vol_correction += (dx_array[j][k] - knick_x[j][k]) * knick_height[j][k] * (bottom_width[j][k] + 0.5 *
					(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
						toe_height_RB[j][k] / tan(toe_angle_RB[j][k])));
			}
			if (bed_tau_c[j][k] == 999 & cohesive_z[j][k] == bed_z[j][k]) {
				skip_bank_erosion_R[j][k] = 2;
				skip_bank_erosion_L[j][k] = 2;
			}
		}
	}

	//update slopes
	double DS_elev = 0;
	for (int j = 0; j < n_nodes; ++j) {
		int index = find(link, j + 1); //Index of DS node
		//Find elevation of downstream XS (either downstream XS in the same reach or the first
		//XS of the next reach)
		for (int k = 0; k < n_xs[j]; ++k) {
			if (k == n_xs[j] - 1) {
				DS_elev = bed_z[index][0];
			}
			else {
				DS_elev = bed_z[j][k + 1];
			}

			//if (k == 0) {
			slope[j][k] = (bed_z[j][k] - DS_elev) / dx_array[j][k];
			//}
			//else {
			//	slope[j][k] = (bed_z[j][k - 1] - DS_elev) / (2.0 * dx_array[j][k]);
			//}
			//Otherwise, check for knickpoints
			if (knick_height[j][k] > 0) {
				slope[j][k] = (bed_z[j][k] - knick_z[j][k]) / knick_x[j][k];
			}
			if (slope[j][k] <= 0) { slope[j][k] = 1e-5; }
		}
	}

	//Print initial outputs
	double time = 0;
	printOutputs(zout_file, bed_z, bottom_width, width_file, time, sinuosity,
		sinuosity_out_file, dx_array, dxout_file, Ds, Ds_log, ps, D50_file, n_xs, stream_power_file,
		slope, q, input_path, slope_file, Rc_file, Rc, top_width);

	//print_bank_outputs(time, bank_sed_mass_file,
	//	bank_p_file, bank_sed_mass, bank_p_mass, cohesive_sed_mass, knick_sed_mass, bed_p_mass,
	//	sed_supply_output_file, sed_supply_output);

	print_geometry(height_LB, toe_height_LB, angle_LB,
		toe_angle_LB, height_RB, toe_height_RB, angle_RB,
		toe_angle_RB, bottom_width, fp_angle, fp_width_R, fp_width_L, bed_z,
		LB_x, geom_file, time);

	print_geometry_final(height_LB, toe_height_LB, angle_LB,
		toe_angle_LB, height_RB, toe_height_RB, angle_RB,
		toe_angle_RB, bottom_width, fp_angle, fp_width_R, fp_width_L, bed_z,
		LB_x, geom_file2, time, n_xs);

	print_knick_outputs(time, knick_out_file, knick_x, dx_array);

	double h_min = 5;
	//Loop through all time steps
	for (long i = 0; i < (n_ts - 1); ++i) {

		Q_iter = (i * dt) / dt_Q;
		//At daily timestep, find unit discharge in channel
		vector<vector<double>> Q_chnl(n_nodes, vector<double>(n_xs_max));
		if (Q_iter != Q_iter_old) {
			for (int j = 0; j < n_nodes; ++j) {
				//Calculate q for each XS
				for (int k = 0; k < n_xs[j]; ++k) {
					//Check discharge and see if it is overbank. Use adjusted discharge if so.
					vector<double> chnl_x(10);
					vector<double> chnl_y(10);
					create_chnl_geom(height_LB[j][k], toe_height_LB[j][k], angle_LB[j][k],
						toe_angle_LB[j][k], height_RB[j][k], toe_height_RB[j][k], angle_RB[j][k],
						toe_angle_RB[j][k], bottom_width[j][k], fp_angle[j][k], fp_width_R[j][k],
						fp_width_L[j][k], bed_z[j][k],
						chnl_x, chnl_y);

					//Check if we shouldn't run bank erosion anymore because toe point equals fp_point
					if (abs(chnl_x[1] - chnl_x[3]) < 1e-5) {
						skip_bank_erosion_L[j][k] = 2;
					}
					if (abs(chnl_x[6] - chnl_x[8]) < 1e-5) {
						skip_bank_erosion_R[j][k] = 2;
					}

					double Q_chnl_val = 0, ws_width = 0;
					typedef pair<double, double> Result;
					Result calc_h = boost::math::tools::brent_find_minima(
						bind(calc_Manning, placeholders::_1, chnl_x, chnl_y, Q[Q_iter][j], slope[j][k], n_chnl[j], n_fp[j],
							Q_chnl_val, ws_width), 0.01, 10.0, 100);
					double h = calc_h.first;

					//Re-call calc_Manning to get Q_chnl and ws_width
					calc_Manning(h, chnl_x, chnl_y, Q[Q_iter][j], slope[j][k], n_chnl[j], n_fp[j], Q_chnl_val,
						ws_width);

					//Don't let calculate Q_chnl exceed total Q (problem for very shallow flows)
					Q_chnl_val = min(Q_chnl_val, Q[Q_iter][j]);

					//if (Q_chnl_val == 0) {
					//	cout << h << " " << slope[j][k] << "\n";
					//}
					//Uses minimum of ws_width and top_width to get water surface width that doesn't extend
					//beyond the channel boundaries
					q[j][k] = Q_chnl_val / min(ws_width, top_width[j][k]);
					Q_chnl[j][k] = Q_chnl_val;

				}
				sed_supply_day[j] = sed_supply[Q_iter][j];
			}
		}

		//Calculate average bottom width
		double DS_width = 0;
		for (int j = 0; j < n_nodes; ++j) {
			int index = find(link, j + 1); //Index of DS node
			for (int k = 0; k < n_xs[j]; ++k) {
				if (k == n_xs[j] - 1) {
					//Only use DS XS width if it is not a tributary junction
					if (index == n_nodes) {
						DS_width = bottom_width[j][k];
					}
					else if (link[1][index] == (n_nodes + 1)) {
						DS_width = bottom_width[index][0];
					}
					else {
						DS_width = bottom_width[j][k];
					}
				}
				else {
					DS_width = bottom_width[j][k + 1];
				}

				avg_width[j][k] = (bottom_width[j][k] + DS_width) / 2.0;
			}
		}

		//Call incision model at each (usually sub-daily) time step
		incision_calcs(n_nodes, i, dt, Q_iter, slope, bed_z, Ds, Ds_log, type, nu, lambda, bottom_width,
			q, link, length, sed_supply_day, dt_output,
			b, height_RB, height_LB, angle_RB, angle_LB, toe_height_RB, toe_height_LB,
			toe_angle_RB, toe_angle_LB, n_xs, dx_array, max_dz, ds_old,
			n_dclass, ps, fs, omegac_star, limit_fun, a, cohesive_z, bed_tau_c, bank_sed, LB_x,
			knick_height, knick_kd, knick_x, knick_z, bank_bedload_prop, bed_sed_out, bed_sed_out_vol,
			fp_angle, fp_width_R, fp_width_L, skip_bank_erosion_L, skip_bank_erosion_R, bank_sed_vol,
			bed_sed_in_vol, dQ_dx_sum, cohesive_bed_vol, transport_factor, cohesive_factor, k_factor,
			cohesive_sed_mass, input_path, sed_supply_output, bed_bedload_prop, cohesive_bedload,
			knick_bedload, avg_width);

		//Call knickpoint model at daily timestep
		if (Q_iter != Q_iter_old & i != 0) {
			//Reset knick_bedload to zero
			for (int j = 0; j < n_nodes; ++j) {
				for (int k = 0; k < n_xs[j]; ++k) {
					knick_bedload[j][k] = 0;
				}
			}

			//Calculate average bottom width
			double DS_width = 0;
			for (int j = 0; j < n_nodes; ++j) {
				int index = find(link, j + 1); //Index of DS node
				for (int k = 0; k < n_xs[j]; ++k) {
					if (k == n_xs[j] - 1) {
						//Only use DS XS width if it is not a tributary junction
						if (index == n_nodes) {
							DS_width = bottom_width[j][k];
						}
						else if (link[1][index] == (n_nodes + 1)) {
							DS_width = bottom_width[index][0];
						}
						else {
							DS_width = bottom_width[j][k];
						}
					}
					else {
						DS_width = bottom_width[j][k + 1];
					}

					avg_width[j][k] = (bottom_width[j][k] + DS_width) / 2.0;
				}
			}

			//Knickpoint model
			KnickErosion(knick_height, knick_kd, knick_x, knick_z, bed_z, Q_chnl, n_nodes, n_xs,
				link, dx_array, height_RB, height_LB, toe_height_RB, toe_height_LB, toe_angle_RB, toe_angle_LB,
				avg_width, knick_bed_vol, knick_sed_mass, cohesive_z, bed_tau_c, bed_bedload_prop,
				knick_bedload, lambda, dt_Q, bed_sed_out_vol, input_path);

			//Calculate p loading from cohesive bed erosion
			for (int j = 0; j < n_nodes; ++j) {
				for (int k = 0; k < n_xs[j]; ++k) {
					bed_p_mass[j][k] = (cohesive_sed_mass[j][k] + knick_sed_mass[j][k]) * bed_p_conc[j] / 1000000.0; //kg P
				}
			}

			time = dt * i / dt_Q;
			print_bank_outputs(time, loading_file, bank_sed_mass, bank_p_mass, cohesive_sed_mass, knick_sed_mass, 
				bed_p_mass);

			//Reset cohesive_sed_mass and sed_supply_output to zero
			for (int j = 0; j < n_nodes; ++j) {
				sed_supply_output[j] = 0;
				for (int k = 0; k < n_xs[j]; ++k) {
					cohesive_sed_mass[j][k] = 0;
				}
			}
		}

		//Only call bank erosion model at daily time step, if requested
		if ((bank_erosion != "none") & (Q_iter != Q_iter_old) & (i != 0)) {
			//if (bank_erosion != "none"){
			BankErosion(tau_c, bank_k, q, bottom_width, top_width, slope, dt_Q, n_nodes,
				i, height_RB, height_LB, angle_RB, angle_LB, toe_height_RB, toe_height_LB, toe_angle_RB,
				toe_angle_LB, cohesion, phi, weight, cohesion_toe, phi_toe, weight_toe, bank_erosion,
				n_xs, bank_sed, dx_array, Rc, sinuosity, n_bends, meandering, bank_bedload_prop,
				Q_iter, Q_iter_old, bank_sed_mass, bank_p_mass, p_conc, LB_x, fp_angle, fp_width_R, fp_width_L,
				bank_tank_RB, bank_tank_LB, skip_bank_erosion_L, skip_bank_erosion_R, lambda, bank_wash_vol,
				fluvial_factor, k_factor, meander_erosion);
		}

		//Output data at user specified frequency
		if ((i != 0) & ((dt * i) % (dt_output * dt_Q) == 0)) {
			time = dt * i / dt_Q;
			printOutputs(zout_file, bed_z, bottom_width, width_file, time, sinuosity,
				sinuosity_out_file, dx_array, dxout_file, Ds, Ds_log, ps, D50_file, n_xs, stream_power_file,
				slope, q, input_path, slope_file, Rc_file, Rc, top_width);

			print_geometry(height_LB, toe_height_LB, angle_LB,
				toe_angle_LB, height_RB, toe_height_RB, angle_RB,
				toe_angle_RB, bottom_width, fp_angle, fp_width_R, fp_width_L, bed_z,
				LB_x, geom_file, time);

			print_geometry_final(height_LB, toe_height_LB, angle_LB,
				toe_angle_LB, height_RB, toe_height_RB, angle_RB,
				toe_angle_RB, bottom_width, fp_angle, fp_width_R, fp_width_L, bed_z,
				LB_x, geom_file2, time, n_xs);

			print_knick_outputs(time, knick_out_file, knick_x, dx_array);
		}

		//Output bed sediment output on daily timestep
		if (MC == 0) {
			if (Q_iter != Q_iter_old) {
				bed_sed_mass_file << bed_sed_out << "\n";
				//bed_sed_mass_file2 << bed_sed_out2 << "\n";
				//bank_sed_vol_file << bank_sed_out << "\n";
				//bed_sed_in_file << bed_sed_in << "\n";
				//dQ_dx_file << dQ_dx_sum << "\n";
				//bank_sed_file << bank_wash_vol << "\n";

				sediment_volume_file << bed_sed_out_vol << "\t" << bed_sed_in_vol << "\t" << "NA" << "\t" <<
					bank_wash_vol << "\t" << cohesive_bed_vol << "\t" << knick_bed_vol << "\t" <<
					bank_sed_vol << "\t" << "NA" << "\t" << dQ_dx_sum << "\n";
			}
			//Output day of simulation if not doing a Monte Carlo simulation

			//if (day != day_old) {
			//	cout << "Day " << day << "\n";
			//}
			if ((i / (n_ts / 20)) != ((i - 1) / (n_ts / 20))) {
				cout << "=";
			}
		}

		//Update day
		Q_iter_old = Q_iter;

	}

	time = dt * (n_ts - 1) / dt_Q;
	print_geometry_final(height_LB, toe_height_LB, angle_LB,
		toe_angle_LB, height_RB, toe_height_RB, angle_RB,
		toe_angle_RB, bottom_width, fp_angle, fp_width_R, fp_width_L, bed_z,
		LB_x, geom_file2, time, n_xs);

	printOutputs(zout_file, bed_z, bottom_width, width_file, time, sinuosity,
		sinuosity_out_file, dx_array, dxout_file, Ds, Ds_log, ps, D50_file, n_xs, stream_power_file,
		slope, q, input_path, slope_file, Rc_file, Rc, top_width);

	print_geometry(height_LB, toe_height_LB, angle_LB,
		toe_angle_LB, height_RB, toe_height_RB, angle_RB,
		toe_angle_RB, bottom_width, fp_angle, fp_width_R, fp_width_L, bed_z,
		LB_x, geom_file, time);

	print_knick_outputs(time, knick_out_file, knick_x, dx_array);

	print_bank_outputs(time, loading_file, bank_sed_mass, bank_p_mass, cohesive_sed_mass, knick_sed_mass,
		bed_p_mass);

	//Subtract eroded knickpoint volume based on remainder of XS
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			if (knick_height[j][k] != 0) {
				knick_vol_correction -= (dx_array[j][k] - knick_x[j][k]) * knick_height[j][k] * (bottom_width[j][k] + 0.5 *
					(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
						toe_height_RB[j][k] / tan(toe_angle_RB[j][k])));
			}
		}
	}

	if (MC == 0){
		bed_sed_mass_file << bed_sed_out << "\n" << flush;

		//Output final bank tank
		double bank_tank_vol = 0;
		for (int j = 0; j < n_nodes; j++) {
			for (int k = 0; k < n_xs[j]; k++) {
				bank_tank_vol += (bank_tank_LB[j][k] + bank_tank_RB[j][k]) * dx_array[j][k];
			}
		}

		sediment_volume_file << bed_sed_out_vol << "\t" << bed_sed_in_vol << "\t" << bank_tank_vol << "\t" <<
			bank_wash_vol << "\t" << cohesive_bed_vol << "\t" << knick_bed_vol << "\t" <<
			bank_sed_vol << "\t" << knick_vol_correction << "\t" << dQ_dx_sum << "\n" << flush;
	}

	//ofstream meander_file(input_path + "Output meander erosion.txt");
	//for (int j = 0; j < n_nodes; ++j) {
	//	for (int k = 0; k < n_xs_max; ++k) {
	//		meander_file << meander_erosion[j][k] << "	";
	//	}
	//	meander_file << "\n";
	//}

	zout_file << flush;
	width_file << flush;
	dxout_file << flush;
	D50_file << flush;
	loading_file << flush;
	geom_file << flush;

	zout_file.close();
	width_file.close();
	dxout_file.close();
	D50_file.close();
	loading_file.close();
	geom_file.close();

	//bank_sed_out_file.close();
	sinuosity_out_file.close();

}

/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//Incision calcs modeling sed transport by size class
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
	vector<vector<double>>& knick_z, vector<double>& bank_bedload_prop, double& bed_sed_out,
	double& bed_sed_out_vol, vector<vector<double>>& fp_angle, vector<vector<double>>& fp_width_R, 
	vector<vector<double>>& fp_width_L, vector<vector<int>>& skip_bank_erosion_L, vector<vector<int>>& skip_bank_erosion_R,
	double& bank_sed_vol, double& bed_sed_in_vol, double& dQ_dx_sum, double& cohesive_bed_vol, double& transport_factor,
	double& cohesive_factor, double& k_factor, vector<vector<double>>& cohesive_sed_mass, string& input_path,
	vector<double>& sed_supply_output, vector<double>& bed_bedload_prop, vector<vector<double>>& cohesive_bedload,
	vector<vector<double>>& knick_bedload, vector<vector<double>>& avg_width) {

	double omega, omega_c;
	int n_xs_max = bed_z[0].size();
	vector<vector<double>> Qb(n_nodes + 1, vector<double>(n_xs_max)); //sediment transport rates (m^3/s)
	const double gamma = 9810; //unit weight of water (N/m^3)
	double delta_z; //change in bed elevation (m)
	double DS_elev; //bed elevation of downstream XS (m)
	double D50 = 0, omegac_star_adj = 0, Qi_US = 0, Qi_DS = 0, ri = 0, ri2 = 0, Psi1 = 0, Psi2 = 0;
	double alpha = 0.5; // weighting for subsurface exchange fraction
	double sg = 2.65; //specific gravity of sediment
	vector<vector<double>> ds(n_nodes, vector<double>(n_xs_max)); //active layer thickness (m)
	vector<vector<double>> Qe_T(n_nodes, vector<double>(n_xs_max)); //total downwind sed transport
	vector<vector<double>> delta_QT(n_nodes, vector<double>(n_xs_max)); //total dQ/dx
	vector<vector<double>> Qsed_out(n_nodes, vector<double>(n_xs_max)); //actual transport out (total)
	vector<vector<double>> Qsed_in_US(n_nodes, vector<double>(n_xs_max)); //actual transport out (total)
	vector<double> fl(n_dclass); //subsurface exchange fraction
	vector<double> ps_val(n_dclass); //vector for storing ps fractions for a XS

	//3D arrays
	typedef boost::multi_array<double, 3 > array;
	array Qi(boost::extents[n_nodes + 1][n_xs_max][n_dclass]); //local transport rate of each size fraction
	array Qe_i(boost::extents[n_nodes + 1][n_xs_max][n_dclass]); //"downwind" transport rate of each size fraction
	array delta_Qi(boost::extents[n_nodes][n_xs_max][n_dclass]); //difference in transport rate between XS for each size fraction

	//Acutal sediment transport out of XS (minimum of capacity and available sediment), by grain size
	array Qsed_outi(boost::extents[n_nodes + 1][n_xs_max][n_dclass]);
	array Qsed_in_USi(boost::extents[n_nodes + 1][n_xs_max][n_dclass]);

	//update slopes
	for (int j = 0; j < n_nodes; ++j) {
		int index = find(link, j + 1); //Index of DS node
		for (int k = 0; k < n_xs[j]; ++k) {
			if (k == n_xs[j] - 1) {
				DS_elev = bed_z[index][0];
			}
			else {
				DS_elev = bed_z[j][k + 1];
			}

			//if (k == 0) {
				slope[j][k] = (bed_z[j][k] - DS_elev) / dx_array[j][k];
			//}
			//else {
			//	slope[j][k] = (bed_z[j][k - 1] - DS_elev) / (2.0 * dx_array[j][k]);
			//}
			//Otherwise, check for knickpoints
			if (knick_height[j][k] > 0) {
				slope[j][k] = (bed_z[j][k] - knick_z[j][k]) / knick_x[j][k];
				//if (j == 4 & k == 0 & knick_x[j][k] < 2) {
				//	cout << "HERE\n";
				//}
			}

			if (slope[j][k] <= 0) { 
				slope[j][k] = 1e-5; 
			}

		}
	}
	
	//Loop through all reaches and XS, calculating sediment transport rate
	for (int j = 0; j < n_nodes; ++j) {
		double US_sed_supply = 0;
		for (int k = 0; k < n_xs[j]; ++k) {
			//Specific stream power
			omega = gamma * slope[j][k] * q_array[j][k];

			for (int m = 0; m < n_dclass; ++m) {
				ps_val[m] = ps[j][k][m];
			}

			//Calculate D50 and active layer thickness (3 * D90)
			D50 = calc_Dx(Ds, ps_val, n_dclass, 0.5, Ds_log, input_path);
			ds[j][k] = min(3.0 * calc_Dx(Ds, ps_val, n_dclass, 0.9, Ds_log, input_path),
				bed_z[j][k] - cohesive_z[j][k]);

			//If first time step, set ds_old to ds
			if (i == 0) { ds_old[j][k] = ds[j][k]; }

			//Loop through all grain sizes
			for (int m = 0; m < n_dclass; ++m) {
				//Adjust sg based on percentage dacite and pumice
				//sg = 2.36 * pow(Ds[m] * 1000, -0.036) * 0.3 + 1.44 * pow(Ds[m] * 1000, -0.09) * 0.3 + 2.65 * 0.4;
				//Find critical stream power for each grain size, adjusting for hiding/protrusion
				omegac_star_adj = omegac_star * pow(Ds[m] / D50, b);
				omega_c = omegac_star_adj * 1000.0 * (9.81 * (sg - 1) * Ds[m]) * sqrt(9.81 * (sg - 1) * Ds[m]);

				//get sediment transport rates at each XS
				Qi[j][k][m] = qs_calc(omega, omega_c, Ds[m], q_array[j][k], ps_val[m], type, 
					transport_factor) * bottom_width[j][k];

				//if (i == (63 * 86400 / dt)) {
				//	cout << Ds[m] << " " << Qi[j][k][m] / bottom_width[j][k] << " " << omega_c << " " << ps_val[m] << "\n";
				//}

				if (isinf(Qi[j][k][m])) {
					cout << "GAAAH!";
				}
				//For most upstream XS, calculate sediment transport capacity if sed_supply = -1
				//if (k == 0 & sed_supply[j] == -1) {
				//	US_sed_supply += qs_calc(omega, omega_c, Ds[m], q_array[j][k], fs[j][k][m], type, transport_factor);
				//}
			}
		}
		//Set sed_supply to US_sed_supply
		//if (sed_supply[j] == -1) { sed_supply[j] = US_sed_supply; }
	}


	//Set last value of Qi to zero for when there are no upstream nodes. This means that there is zero
	//incoming sediment from the upstream channel because there is no upstream channel
	for (int m = 0; m < n_dclass; ++m) {
		Qi[n_nodes][0][m] = 0;
	}

	//Again, loop through all reaches and XS to calculate the "downwind" adjusted sediment transport rate
	//using a flux limiting function
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			for (int m = 0; m < n_dclass; ++m) {
				if (k == 0) {
					//Most upstream XS in reach
					if (sed_supply[j] < 0) {
						Qi_US = Qi[j][k][m] * abs(sed_supply[j]);
					}
					else {
						Qi_US = sed_supply[j] * fs[j][k][m] * bottom_width[j][k] + Qi[link[0][j] - 1][n_xs[link[0][j] - 1] - 1][m] +
							Qi[link[1][j] - 1][n_xs[link[1][j] - 1] - 1][m];
					}
					//If there isn't any incoming sediment, set transport rate US equal to transport rate at XS
					if (Qi_US == 0) {
						Qi_US = Qi[j][k][m];
					}
					Qi_DS = Qi[j][k + 1][m];
				}
				else if (k == n_xs[j] - 1) {
					//Most downstream XS in reach
					Qi_US = Qi[j][k - 1][m];
					int index = find(link, j + 1); //Index of DS node
					Qi_DS = Qi[index][0][m];
				}
				else {
					//All other XS
					Qi_US = Qi[j][k - 1][m];
					Qi_DS = Qi[j][k + 1][m];
				}

				//Catch the case where the fluxes are equal
				if (Qi_DS == Qi[j][k][m] || Qi[j][k][m] == Qi_US) {
					Qe_i[j][k][m] = Qi[j][k][m];
				}
				else {
					//Use limiter function
					ri = (Qi_DS - Qi[j][k][m]) / (Qi[j][k][m] - Qi_US);
					ri2 = 1.0 / ri;
					Psi1 = Limiter(ri, limit_fun);
					Psi2 = Limiter(ri2, limit_fun);
					Qe_i[j][k][m] = Qi[j][k][m] + 1.0 / 4.0 * ((1 - a) * Psi1 + (1 + a) * ri * Psi2) *
						(Qi[j][k][m] - Qi_US);
				}

				//Don't let transport rate be negative
				if (Qe_i[j][k][m] < 0) { Qe_i[j][k][m] = 0; }

				//Sum total transport rate (sum of all size-fraction transport rates)
				Qe_T[j][k] += Qe_i[j][k][m];

				if (Qe_T[j][k] != Qe_T[j][k]) {
					cout << "GAAAH!";
				}
			}
		}
	}

	//Set last value of Qe_i to zero for when there are no upstream nodes. This means that there is zero
	//incoming sediment from the upstream channel because there is no upstream channel
	for (int m = 0; m < n_dclass; ++m) {
		Qe_i[n_nodes][0][m] = 0;
	}

	//Bed evolution loop
	double QT_US = 0;
	double bed_supply_i = 0;
	double bed_supply = 0;
	//reset Qsed_out to zero
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			Qsed_out[j][k] = 0;
		}
	}

	//Loop through each reach and XS, calculating dQ/dx, change in bed elevation, and updating
	//grain size distribution
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			if (k == 0) {
				//Most upstream XS in segment
				QT_US = 0;
				for (int m = 0; m < n_dclass; ++m) {
					//Calculate available supply of sediment on bed (volume bed of that size fraction,
					//divided by dt to get a rate)
					bed_supply_i  = (((bed_z[j][k] - ds[j][k] - cohesive_z[j][k]) * fs[j][k][m] + ds[j][k] * ps[j][k][m]) *
						avg_width[j][k] * dx_array[j][k] * (1 - lambda) + ((bed_z[j][k] - cohesive_z[j][k]) *
							0.5 * (toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) + toe_height_RB[j][k] /
								tan(toe_angle_RB[j][k])) * fs[j][k][m] * dx_array[j][k] * (1 - lambda))) / dt;

					if (bed_supply_i < 0 | (bed_z[j][k] - cohesive_z[j][k]) == 0) {
						bed_supply_i = 0;
					}

					//Actual amount of sediment coming in (sum of watershed sediment supply (if any),
					//sediment from upstream reach(s), and any sediment from local bank erosion)
					if (sed_supply[j] < 0) {
						Qi_US = Qe_i[j][k][m] * abs(sed_supply[j]) + Qsed_outi[link[0][j] - 1][n_xs[link[0][j] - 1] - 1][m] +
							Qsed_outi[link[1][j] - 1][n_xs[link[1][j] - 1] - 1][m] +
							(bank_sed[j][k] + cohesive_bedload[j][k] + knick_bedload[j][k]) * fs[j][k][m];
						bed_sed_in_vol += Qe_i[j][k][m] * abs(sed_supply[j]) * dt;
						sed_supply_output[j] += Qe_i[j][k][m] * abs(sed_supply[j]) * dt / bottom_width[j][k];
					}
					else {
						Qi_US = sed_supply[j] * fs[j][k][m] * bottom_width[j][k] + Qsed_outi[link[0][j] - 1][n_xs[link[0][j] - 1] - 1][m] +
							Qsed_outi[link[1][j] - 1][n_xs[link[1][j] - 1] - 1][m] +
							(bank_sed[j][k] + cohesive_bedload[j][k] + knick_bedload[j][k]) * fs[j][k][m];
						bed_sed_in_vol += sed_supply[j] * bottom_width[j][k] * dt;
						sed_supply_output[j] += sed_supply[j] * dt * fs[j][k][m];
					}

					//Calculate dQs/dx, using (capacity - sed_in) for each size fraction
					delta_Qi[j][k][m] = (Qe_i[j][k][m] - Qi_US) / dx_array[j][k];
					QT_US += Qi_US;
					Qsed_in_USi[j][k][m] = Qi_US;

					//Sediment out of XS is equal to minimum of capacity and 
					//available bed sediment (bed_supply_i) + sediment coming from upstream
					Qsed_outi[j][k][m] = min(Qe_i[j][k][m], bed_supply_i + Qi_US);
					//Qsed_outi[j][k][m] = min(Qe_i[j][k][m], bed_supply_i);

					//if (bed_supply_i == 0) {
					//	cout << Qsed_outi[j][k][m] << " " << Qe_i[j][k][m] << " " << Qi_US << "\n";
					//}

					//Get cumulative sediment out (sum of all size fractions)
					Qsed_out[j][k] += Qsed_outi[j][k][m];
				}
				Qsed_in_US[j][k] = QT_US;
				//Calculate total dQs/dx
				delta_QT[j][k] = (Qe_T[j][k] - QT_US) / dx_array[j][k];

				}
			else {
				//All other XS
				//Calculate total dQs/dx
				delta_QT[j][k] = (Qe_T[j][k] - Qsed_out[j][k - 1] - bank_sed[j][k] - cohesive_bedload[j][k] - knick_bedload[j][k]) / dx_array[j][k];
				for (int m = 0; m < n_dclass; ++m) {
					//Calculate dQs/dx by size fraction
					delta_Qi[j][k][m] = (Qe_i[j][k][m] - Qsed_outi[j][k - 1][m] - (bank_sed[j][k] + 
						cohesive_bedload[j][k] + knick_bedload[j][k]) * fs[j][k][m]) / dx_array[j][k];

					//Calculate available supply of sediment on bed (volume bed of that size fraction,
					//divided by dt to get a rate)
					bed_supply_i = (((bed_z[j][k] - ds[j][k] - cohesive_z[j][k]) * fs[j][k][m] + ds[j][k] * ps[j][k][m]) *
						avg_width[j][k] * dx_array[j][k] * (1 - lambda) + ((bed_z[j][k] - cohesive_z[j][k]) *
							0.5 * (toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) + toe_height_RB[j][k] / 
								tan(toe_angle_RB[j][k])) * fs[j][k][m] * dx_array[j][k] * (1 - lambda))) / dt;
					
					if (bed_supply_i < 0 | (bed_z[j][k] - cohesive_z[j][k]) == 0) {
						bed_supply_i = 0;
					}
					//Sediment out of XS is equal to minimum of capacity and 
					//available bed sediment (bed_supply_i) + sediment coming from upstream and bank erosion
					Qsed_outi[j][k][m] = min(Qe_i[j][k][m], bed_supply_i + Qsed_outi[j][k - 1][m] + 
						(bank_sed[j][k] + cohesive_bedload[j][k] + knick_bedload[j][k]) * fs[j][k][m]);

					//if (j == 1 & k == 9 & Qsed_outi[j][k][m] != 0 & (Qe_i[j][k][m] != Qsed_outi[j][k - 1][m])) {
					//	cout << Qsed_outi[j][k][m] << " " << Qe_i[j][k][m] << " " << Qsed_outi[j][k - 1][m] << " " << bed_supply_i << "\n";
					//}

					Qsed_out[j][k] += Qsed_outi[j][k][m];
				}
			}

			if (bottom_width[j][k] < 0) {
				cout << "CRAP - width\n";
			}

			//Export total bed material load leaving the watershed
			if ((j == (n_nodes - 1)) & k == (n_xs[j] - 1)) {
				for (int m = 0; m < n_dclass; ++m) {
					//Adjust sg based on percentage dacite and pumice
					//sg = 2.36 * pow(Ds[m] * 1000, -0.036) * 0.3 + 1.44 * pow(Ds[m] * 1000, -0.09) * 0.3 + 2.65 * 0.4;

					bed_sed_out += Qsed_outi[j][k][m] * dt * 1000.0 * sg;
				}
				bed_sed_out_vol += Qsed_out[j][k] * dt;
				
			}
			bank_sed_vol += bank_sed[j][k] * dt;
			
			//EXNER EQUATION FOR CHANGE IN BED ELEVATION
			delta_z = -1 / (1 - lambda) * (delta_QT[j][k]) * dt / avg_width[j][k];

			if (delta_z < 0) {
				double delta_z2 = delta_z * avg_width[j][k] / (avg_width[j][k] + 0.5 *
					(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) + toe_height_RB[j][k] / tan(toe_angle_RB[j][k])));
				delta_z = delta_z2;
			}else if (delta_z > 0){
				vector<double> chnl_x(10);
				vector<double> chnl_y(10);
				double bed_z2 = 0;
				create_chnl_geom(height_LB[j][k], toe_height_LB[j][k], angle_LB[j][k],
					toe_angle_LB[j][k], height_RB[j][k], toe_height_RB[j][k], angle_RB[j][k],
					toe_angle_RB[j][k], bottom_width[j][k], fp_angle[j][k], fp_width_R[j][k],
					fp_width_L[j][k], bed_z2,
					chnl_x, chnl_y);

				typedef pair<double, double> Result;
				Result dz_calc = boost::math::tools::brent_find_minima(
					bind(calc_dz_adj, placeholders::_1, toe_height_RB[j][k], toe_height_LB[j][k], 
						toe_angle_RB[j][k], toe_angle_LB[j][k], angle_RB[j][k], angle_LB[j][k], 
						avg_width[j][k], bank_bedload_prop[j], fp_angle[j][k], fp_width_R[j][k],
						fp_width_L[j][k], height_RB[j][k], height_LB[j][k], bed_z[j][k], delta_z,
						chnl_x, chnl_y), 0.0, delta_z, 100);

				delta_z = dz_calc.first;
				double area_fp, area_chnl, perim, perim_fp, perim_chnl, ws_width;
				double area_calc = calc_flow_area(chnl_x, chnl_y, delta_z, area_fp,
					area_chnl, perim, perim_fp, perim_chnl, ws_width);
				dQ_dx_sum += abs(area_calc + 1 / (1 - lambda) * (delta_QT[j][k]) * dt);
			}

			if (bed_z[j][k] < cohesive_z[j][k]) {
				cout << "HELP!\n";
			}
			/////////////////////////////////////////////////////////////////////////
			//Check if bed is encountering cohesive layer
			double dz_cohesive = 0, dz_alluvium = 0;
			if ((delta_z <= 0) & ((bed_z[j][k] + delta_z <= cohesive_z[j][k]) | (bed_z[j][k] <= cohesive_z[j][k]))) {
				//(re)calculate stream power
				omega = 9810.0 * slope[j][k] * q_array[j][k];

				//check if there is still some alluvium in channel or if it is all cohesive
				if (bed_z[j][k] > cohesive_z[j][k]) {
					dz_alluvium = - (bed_z[j][k] - cohesive_z[j][k]);
				}
				if (dz_alluvium > 0) {
					cout << dz_alluvium << "\n";
				}

				//Call cohesive erosion function
				dz_cohesive = cohesive_incision(omega, bed_tau_c[j][k], dt, cohesive_factor, k_factor);

				//Add two together to get total dz
				delta_z = dz_alluvium + dz_cohesive;

				//Adjust cohesive if non-zero
				//Adjust cohesive bed elevation and calculate volume of eroded cohesive material
				cohesive_z[j][k] += dz_cohesive;
				cohesive_bed_vol += (abs(dz_cohesive) * (avg_width[j][k] + 0.5 *
					(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) + toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))
					* dx_array[j][k]) * (1 - bed_bedload_prop[j]);
				//Cumulative sed mass over the entire day. Then value will be reset to zero for the next day.
				//Must incorporate the full area of bed erosion (including the "adjusted" area from the bank toes)
				cohesive_sed_mass[j][k] += abs(dz_cohesive) * (avg_width[j][k] + 0.5 *
					(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) + toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))
					* dx_array[j][k] * 1600 * (1 - bed_bedload_prop[j]); //kg; assumes bulk soil density of 1.6 Mg/m^3
				//Bedload material from cohesive bed erosion (convert to rate m^3/s)
				cohesive_bedload[j][k] = (abs(dz_cohesive) * (avg_width[j][k] + 0.5 *
					(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) + toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))
					* dx_array[j][k]) * bed_bedload_prop[j] * (1 - lambda) / dt;
			
			}
			else {
				cohesive_bedload[j][k] = 0;
			}
			///////////////////////////////////////////////////////////////////////

			//Update bed elevation
			bed_z[j][k] = bed_z[j][k] + delta_z;

			if (bed_z[j][k] != bed_z[j][k]) {
				cout << "ERROR\n";
			}

			max_dz = max(max_dz, abs(delta_z));

			//Update bank height but make sure it doesn't go negative (i.e. the channel completely fills in)
			double old_h_RB = height_RB[j][k];
			double old_h_LB = height_LB[j][k];
			height_RB[j][k] = height_RB[j][k] - delta_z;
			height_LB[j][k] = height_LB[j][k] - delta_z;

			double RB_h_filled = 0;
			double LB_h_filled = 0;
			if (height_RB[j][k] < 0) {
				RB_h_filled = 1;
				bottom_width[j][k] += (old_h_RB - toe_height_RB[j][k]) / tan(angle_RB[j][k]);
				//Add width from fp
				if (fp_angle[j][k] != 0 & fp_width_R[j][k] != 0 & (abs(height_RB[j][k]) / tan(fp_angle[j][k]) < fp_width_R[j][k])) {
					bottom_width[j][k] += abs(height_RB[j][k]) / tan(fp_angle[j][k]);
					fp_width_R[j][k] -= abs(height_RB[j][k]) / tan(fp_angle[j][k]);
					height_RB[j][k] = fp_width_R[j][k] * tan(fp_angle[j][k]);
					angle_RB[j][k] = fp_angle[j][k];
					fp_width_R[j][k] = 0;
					skip_bank_erosion_R[j][k] = 1;
				}
				else {
					bottom_width[j][k] += fp_width_R[j][k];
					fp_width_R[j][k] = 0;
					height_RB[j][k] = 10;
					angle_RB[j][k] = M_PI / 2.0;
					skip_bank_erosion_R[j][k] = 2;
				}				
			}
			if (height_LB[j][k] < 0) {
				LB_h_filled = 1;
				bottom_width[j][k] += (old_h_LB - toe_height_LB[j][k]) / tan(angle_LB[j][k]);
				LB_x[j][k] -= (old_h_LB - toe_height_LB[j][k]) / tan(angle_LB[j][k]);
				//Add width from fp
				if (fp_angle[j][k] != 0 & fp_width_L[j][k] != 0 & (abs(height_LB[j][k]) / tan(fp_angle[j][k]) < fp_width_L[j][k])) {
					bottom_width[j][k] += abs(height_LB[j][k]) / tan(fp_angle[j][k]);
					LB_x[j][k] -= abs(height_LB[j][k]) / tan(fp_angle[j][k]);
					fp_width_L[j][k] -= abs(height_LB[j][k]) / tan(fp_angle[j][k]);
					height_LB[j][k] = fp_width_L[j][k] * tan(fp_angle[j][k]);
					angle_LB[j][k] = fp_angle[j][k];
					fp_width_L[j][k] = 0;
					skip_bank_erosion_L[j][k] = 1;
				}
				else {
					bottom_width[j][k] += fp_width_L[j][k];
					LB_x[j][k] -= fp_width_L[j][k];
					fp_width_L[j][k] = 0;
					height_LB[j][k] = 10;
					angle_LB[j][k] = M_PI / 2.0;
					skip_bank_erosion_L[j][k] = 2;
				}
			}

			//Update toe height the same way
			toe_height_RB[j][k] = toe_height_RB[j][k] - delta_z;
			toe_height_LB[j][k] = toe_height_LB[j][k] - delta_z;

			double toe_angle_RB_old = toe_angle_RB[j][k];
			double toe_angle_LB_old = toe_angle_LB[j][k];
			//Update toe angle - ONLY IF INCISING IF WIDTH CHANGES WHILE AGGRADING
			//if (delta_z < 0) {
				toe_angle_RB[j][k] = atan(toe_height_RB[j][k] / (toe_height_RB[j][k] + delta_z) *
					tan(toe_angle_RB[j][k]));
				toe_angle_LB[j][k] = atan(toe_height_LB[j][k] / (toe_height_LB[j][k] + delta_z) *
					tan(toe_angle_LB[j][k]));
			//}

			//If toe height is less than zero (it filled in) then make toe height half of bank height
			//and set their angles equal
			if (toe_height_RB[j][k] < 0 & toe_height_LB[j][k] < 0) {
				//both toe's filled in
				bottom_width[j][k] += (toe_height_RB[j][k] + delta_z) / tan(toe_angle_RB_old);
				if (RB_h_filled == 0) {
					bottom_width[j][k] += abs(toe_height_RB[j][k]) / tan(angle_RB[j][k]);
				}
				toe_height_RB[j][k] = height_RB[j][k] / 2.0;
				toe_angle_RB[j][k] = angle_RB[j][k];
				bottom_width[j][k] += (toe_height_LB[j][k] + delta_z) / tan(toe_angle_LB_old);
				LB_x[j][k] -= (toe_height_LB[j][k] + delta_z) / tan(toe_angle_LB_old);
				if (LB_h_filled == 0) {
					bottom_width[j][k] += abs(toe_height_LB[j][k]) / tan(angle_LB[j][k]);
					LB_x[j][k] -= abs(toe_height_LB[j][k]) / tan(angle_LB[j][k]);
				}
				toe_height_LB[j][k] = height_LB[j][k] / 2.0;
				toe_angle_LB[j][k] = angle_LB[j][k];
			}else if (toe_height_RB[j][k] < 0){
				//only RB filled in
				bottom_width[j][k] += (toe_height_RB[j][k] + delta_z) / tan(toe_angle_RB_old);
				if (RB_h_filled == 0) {
					bottom_width[j][k] += abs(toe_height_RB[j][k]) / tan(angle_RB[j][k]);
				}
				toe_height_RB[j][k] = height_RB[j][k] / 2.0;
				toe_angle_RB[j][k] = angle_RB[j][k];

				//Widen channel based on LB - INCLUDE IF WIDTH CHANGES WITH AGGRADATION
				//bottom_width[j][k] += delta_z / tan(toe_angle_LB[j][k]);
				//LB_x[j][k] -= delta_z / tan(toe_angle_LB[j][k]);
			}else if (toe_height_LB[j][k] < 0) {
				//only LB filled in
				bottom_width[j][k] += (toe_height_LB[j][k] + delta_z) / tan(toe_angle_LB_old);
				LB_x[j][k] -= (toe_height_LB[j][k] + delta_z) / tan(toe_angle_LB_old);
				if (LB_h_filled == 0) {
					bottom_width[j][k] += abs(toe_height_LB[j][k]) / tan(angle_LB[j][k]);
					LB_x[j][k] -= abs(toe_height_LB[j][k]) / tan(angle_LB[j][k]);
				}
				toe_height_LB[j][k] = height_LB[j][k] / 2.0;
				toe_angle_LB[j][k] = angle_LB[j][k];

				//Widen channel based on RB- INCLUDE IF WIDTH CHANGES WITH AGGRADATION
				//bottom_width[j][k] += delta_z / tan(toe_angle_RB[j][k]);
			}
			/*else if (delta_z > 0) { //- INCLUDE IF WIDTH CHANGES WITH AGGRADATION
				//Neither toe is filled in, but still aggrading
				bottom_width[j][k] += delta_z * (1 / tan(toe_angle_RB[j][k]) + 1 / tan(toe_angle_LB[j][k]));

				LB_x[j][k] -= delta_z / tan(toe_angle_LB[j][k]);

				//if (j == 15 & k == 0) {
				//	cout << "HELP!";
				//}
			}*/
			///////////////////////////////////////////////////////////////////////////
			
			if (bottom_width[j][k] < 0) {
				cout << "CRAP - width\n";
			}

			//////////////////////////////////////////////////////////////////////////
			//Update grain size fractions - if ds is not zero
			//if (ds[j][k] > 0) {
				//if (j == 0 & k == 4) {
				//	cout << "HEFE\n";
				//}
				double ps_sum = 0;
				for (int m = 0; m < n_dclass; ++m) {
					if (delta_z < 0) {
						//interface exchange fraction equals the subsurface portion if degrading
						fl[m] = fs[j][k][m];
					}
					else {
						//interface exchange fraction is a weighted average of the bedload and bed material
						//if the channel is aggrading
						if (Qsed_out[j][k] == 0) {
							fl[m] = ps[j][k][m];
						}
						else {
							fl[m] = alpha * ps[j][k][m] + (1 - alpha) * Qsed_outi[j][k][m] / Qsed_out[j][k];
						}
					}

					//Get grain size fraction so we don't have to keep referencing the large 3-D array
					double ps_a = ps[j][k][m];

					//Update grain size fraction
					ps_a = ps_a + 1 / ds[j][k] * (fl[m] - ps_a) * (ds_old[j][k] -
						ds[j][k]) + 1 / (ds[j][k] * (1 - lambda) * avg_width[j][k]) *
						(-delta_Qi[j][k][m] + fl[m] * delta_QT[j][k]) * dt;

					//If aggrading onto a cohesive bed, set new grain size fraction equal to the bedload fraction
					if (ds[j][k] == 0) {
						if (k != 0) {
							if (Qsed_out[j][k - 1] != 0) {
								ps_a = Qsed_outi[j][k - 1][m] / Qsed_out[j][k - 1];
							}
							else {
								ps_a = ps[j][k][m];
							}
						}
						else {
							if (Qsed_in_US[j][k] != 0){
								ps_a = Qsed_in_USi[j][k][m] / Qsed_in_US[j][k];
							}
							else {
								ps_a = ps[j][k][m];
							}
						}
					}

					//Is ps negative? Set threshold to a really low value (1e-5) because of rounding errors
					if (ps_a < 1e-5) { ps_a = 0.0; }

					//Sum new grain size fractions
					ps_sum += ps_a;

					//Update value in 3-D array
					ps[j][k][m] = ps_a;
				}

				//Exit with error if ps_sum = 0
				if (ps_sum == 0) { exit(1); }

				//Renormalize to 0-1 scale
				for (int m = 0; m < n_dclass; ++m) {
					ps[j][k][m] /= ps_sum;
				}
			//}
			//Set ds_old to ds (active layer thickness)
			ds_old[j][k] = ds[j][k];

			//////////////////////////////////////////////////////////////////////////////////////////
		}
	}
}

double cohesive_incision(double& omega, double& bed_tau_c, int& dt, double& cohesive_factor, double& k_factor) {
	double E = 0; //eroded distance (m)

	//calculate bed shear stress from stream power
	double tau = 1.96 * pow(omega, 0.72) * cohesive_factor;

	//empirical relationship to get K from tau_c (from Simon et al., 2010)
	double bed_k = 1.6 * pow(bed_tau_c, -0.826) * 1e-6 * k_factor;

	//If shear stress exceeds the critical value, calculate erosion, otherwise set erosion to zero
	if (tau > bed_tau_c) {
		E = - bed_k * dt * (tau - bed_tau_c);
	}
	else {
		E = 0.0;
	}

	return(E);

}

void KnickErosion(vector<vector<double>>& knick_height, vector<vector<double>>& knick_kd, vector<vector<double>>& knick_x,
	vector<vector<double>>& knick_z, vector<vector<double>>& bed_z, vector<vector<double>>& Q_chnl, int& n_nodes, 
	vector<int>& n_xs, vector<vector<int>>& link, vector<vector<double>>& dx_array, vector<vector<double>>& height_RB, 
	vector<vector<double>>& height_LB, vector<vector<double>>& toe_height_RB, vector<vector<double>>& toe_height_LB,
	vector<vector<double>>& toe_angle_RB, vector<vector<double>>& toe_angle_LB, 
	vector<vector<double>>& avg_width, double& knick_bed_vol, vector<vector<double>>& knick_sed_mass,
	vector<vector<double>>& cohesive_z, vector<vector<double>>& bed_tau_c, vector<double>& bed_bedload_prop,
	vector<vector<double>>& knick_bedload, double& lambda, int& dt, double& bed_sed_out_vol, string& input_path) {
	
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			//Reset knick sed mass as zero
			knick_sed_mass[j][k] = 0;
			double knick_vol = 0;
			//Only do anything if there is a knickpoing
			if (knick_height[j][k] != 0) {
				//Get daily cumulative flow (m^3)
				double Q = Q_chnl[j][k] * dt; //m^3/s * s/day
				//Calculate Ehc from kd
				//double Ehc = 17.8 + (16.5 * knick_kd[j][k]); //No vegetation effects - results in really fast migration
				double Ehc = 17.8 + (16.5 * knick_kd[j][k] - 15 * 1.4);
				if (Ehc < 0) {
					ofstream error_file(input_path + "MODEL ERRORS.txt");
					error_file << "Negative value for knickpoint erosion. Increase Kd (minimum acceptable value is 3.2/16.5 or about 0.194.\n";
					exit(1);
				}
				//Knickpoint migration distance
				double L = 0.00126 * Ehc * sqrt(Q) * pow(knick_height[j][k], 0.225);
				if (L == 0) {
					cout << L << " " << Q << " " << Ehc << " " << knick_height[j][k] << "\n";
				}
				//Update variables
				if (L > knick_x[j][k]) {

					if (L - knick_x[j][k] > 100) {
						cout << "HEY!\n";
					}
					//knickpoint migrated past US XS
					if (k == 0 & (link[0][j] == link[1][j])) {
					//	//headwater tributary - update bed_z of XS and get rid of knickpoint (done after if statements

					}
					else if (k == 0) {
						//Non-headwater reach - send knickpoint up both tributaries
						int trib1 = link[0][j] - 1;
						int trib2 = link[1][j] - 1;

						//Trib 1
						if (trib1 != n_nodes) {
							if (bed_tau_c[trib1][n_xs[trib1] - 1] != 999){
								knick_height[trib1][n_xs[trib1] - 1] = knick_height[j][k];
								knick_x[trib1][n_xs[trib1] - 1] = dx_array[trib1][n_xs[trib1] - 1] - L - knick_x[j][k];
								knick_z[trib1][n_xs[trib1] - 1] = bed_z[trib1][n_xs[trib1] - 1] - 
									(bed_z[trib1][n_xs[trib1] - 1] - bed_z[j][k]) /
									dx_array[trib1][n_xs[trib1] - 1] * knick_x[trib1][n_xs[trib1] - 1];
								knick_kd[trib1][n_xs[trib1] - 1] = knick_kd[j][k];
								knick_vol += ((L - knick_x[j][k]) * knick_height[j][k] *
									(avg_width[trib1][n_xs[trib1] - 1] + 0.5 *
									(toe_height_LB[trib1][n_xs[trib1] - 1] / tan(toe_angle_LB[trib1][n_xs[trib1] - 1]) +
										toe_height_RB[trib1][n_xs[trib1] - 1] / tan(toe_angle_RB[trib1][n_xs[trib1] - 1])))) *
									(1 - bed_bedload_prop[trib1]);
								knick_bedload[j][k] += ((L - knick_x[j][k]) * knick_height[j][k] *
									(avg_width[trib1][n_xs[trib1] - 1] + 0.5 *
									(toe_height_LB[trib1][n_xs[trib1] - 1] / tan(toe_angle_LB[trib1][n_xs[trib1] - 1]) +
										toe_height_RB[trib1][n_xs[trib1] - 1] / tan(toe_angle_RB[trib1][n_xs[trib1] - 1])))) *
										(bed_bedload_prop[trib1]) * (1 - lambda) / dt;
							}
						}

						//Trib 2
						if (trib2 != n_nodes) {
							if (bed_tau_c[trib2][n_xs[trib2] - 1] != 999) {
								knick_height[trib2][n_xs[trib2] - 1] = knick_height[j][k];
								knick_x[trib2][n_xs[trib2] - 1] = dx_array[trib2][n_xs[trib2] - 1] - L - knick_x[j][k];
								knick_z[trib2][n_xs[trib2] - 1] = bed_z[trib2][n_xs[trib2] - 1] - 
									(bed_z[trib2][n_xs[trib2] - 1] - bed_z[j][k]) /
									dx_array[trib2][n_xs[trib2] - 1] * knick_x[trib2][n_xs[trib2] - 1];
								knick_kd[trib2][n_xs[trib2] - 1] = knick_kd[j][k];
								knick_vol += ((L - knick_x[j][k]) * knick_height[j][k] *
									(avg_width[trib2][n_xs[trib2] - 1] + 0.5 *
									(toe_height_LB[trib2][n_xs[trib2] - 1] / tan(toe_angle_LB[trib2][n_xs[trib2] - 1]) +
										toe_height_RB[trib2][n_xs[trib2] - 1] / tan(toe_angle_RB[trib2][n_xs[trib2] - 1])))) *
										(1 - bed_bedload_prop[trib2]);
								knick_bedload[j][k] += ((L - knick_x[j][k]) * knick_height[j][k] *
									(avg_width[trib2][n_xs[trib2] - 1] + 0.5 *
									(toe_height_LB[trib2][n_xs[trib2] - 1] / tan(toe_angle_LB[trib2][n_xs[trib2] - 1]) +
										toe_height_RB[trib2][n_xs[trib2] - 1] / tan(toe_angle_RB[trib2][n_xs[trib2] - 1])))) *
										(bed_bedload_prop[trib2]) * (1 - lambda) / dt;
							}
						}
					}
					else if (bed_tau_c[j][k - 1] != 999) {
						//Knickpoint stays in same reach - doesn't migrate US if US XS has grade control (tau_c = 999)
						knick_height[j][k - 1] = knick_height[j][k];
						knick_x[j][k - 1] = dx_array[j][k - 1] - L - knick_x[j][k];
						knick_z[j][k - 1] = bed_z[j][k - 1] - (bed_z[j][k - 1] - bed_z[j][k]) / dx_array[j][k - 1] * knick_x[j][k - 1];
						knick_kd[j][k - 1] = knick_kd[j][k];
						knick_vol += ((L - knick_x[j][k]) * knick_height[j][k] * (avg_width[j][k - 1] + 0.5 *
							(toe_height_LB[j][k - 1] / tan(toe_angle_LB[j][k - 1]) +
								toe_height_RB[j][k - 1] / tan(toe_angle_RB[j][k - 1])))) * (1 - bed_bedload_prop[j]);
						knick_bedload[j][k] += ((L - knick_x[j][k]) * knick_height[j][k] * (avg_width[j][k - 1] + 0.5 *
							(toe_height_LB[j][k - 1] / tan(toe_angle_LB[j][k - 1]) +
								toe_height_RB[j][k - 1] / tan(toe_angle_RB[j][k - 1])))) * (bed_bedload_prop[j]) * (1 - lambda) / dt;
					}

					//Update bank heights/angles, bed elevations, and set knickpoint variables to zero
					knick_vol += (knick_x[j][k] * knick_height[j][k] * (avg_width[j][k] + 0.5 *
						(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
							toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (1 - bed_bedload_prop[j]);
					if (k == n_xs[j] - 1 & j == n_nodes - 1) {
						bed_sed_out_vol += (knick_x[j][k] * knick_height[j][k] * (avg_width[j][k] + 0.5 *
							(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
								toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (bed_bedload_prop[j]) * (1 - lambda);
					}
					else if (j != n_nodes - 1) {
						int index = find(link, j + 1); //Index of DS node
						knick_bedload[index][0] += (knick_x[j][k] * knick_height[j][k] * (avg_width[j][k] + 0.5 *
							(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
								toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (bed_bedload_prop[j]) * (1 - lambda) / dt;
					}
					else{
						knick_bedload[j][k + 1] += (knick_x[j][k] * knick_height[j][k] * (avg_width[j][k] + 0.5 *
							(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
								toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (bed_bedload_prop[j]) * (1 - lambda) / dt;
					}
					height_RB[j][k] += knick_height[j][k];
					height_LB[j][k] += knick_height[j][k];
					toe_height_RB[j][k] += knick_height[j][k];
					toe_height_LB[j][k] += knick_height[j][k];
					toe_angle_RB[j][k] = atan(toe_height_RB[j][k] * tan(toe_angle_RB[j][k]) /
						(toe_height_RB[j][k] - knick_height[j][k]));
					toe_angle_LB[j][k] = atan(toe_height_LB[j][k] * tan(toe_angle_LB[j][k]) /
						(toe_height_LB[j][k] - knick_height[j][k]));
					bed_z[j][k] -= knick_height[j][k];
					if (cohesive_z[j][k] > bed_z[j][k]) { cohesive_z[j][k] = bed_z[j][k]; }
					//cohesive_z[j][k] = bed_z[j][k];
					knick_height[j][k] = 0;
					knick_x[j][k] = 0;
					knick_z[j][k] = 0;
					knick_kd[j][k] = 0;
				}
				else {
					//Knickpoint just moved within same XS - just update elevation and x
					knick_x[j][k] -= L;
					knick_z[j][k] = bed_z[j][k] - (bed_z[j][k] - knick_z[j][k]) / (knick_x[j][k] + L) * knick_x[j][k];
					knick_vol += (L * knick_height[j][k] * (avg_width[j][k] + 0.5 *
						(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
							toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (1 - bed_bedload_prop[j]);
					if (k == n_xs[j] - 1 & j == n_nodes - 1) {
						bed_sed_out_vol += (L * knick_height[j][k] * (avg_width[j][k] + 0.5 *
							(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
								toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (bed_bedload_prop[j]) * (1 - lambda);
					}
					else if (j != n_nodes - 1) {
						int index = find(link, j + 1); //Index of DS node
						knick_bedload[index][0] += (L * knick_height[j][k] * (avg_width[j][k] + 0.5 *
							(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
								toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (bed_bedload_prop[j]) * (1 - lambda) / dt;
					}
					else {
						knick_bedload[j][k + 1] += (L * knick_height[j][k] * (avg_width[j][k] + 0.5 *
							(toe_height_LB[j][k] / tan(toe_angle_LB[j][k]) +
								toe_height_RB[j][k] / tan(toe_angle_RB[j][k])))) * (bed_bedload_prop[j]) * (1 - lambda) / dt;
					}
				}
				knick_bed_vol += knick_vol;
				knick_sed_mass[j][k] = knick_vol * 1600; //kg; assumes soil bulk density of 1.6 Mg/kg

				if (knick_vol <= 0) {
					cout << "AHEH\n";
				}
			}
		}
	}

}