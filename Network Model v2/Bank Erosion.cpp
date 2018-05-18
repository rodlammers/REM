
#include "stdafx.h"
#include "Bank Erosion.h"
//using namespace std;

//BankShear calculates eroded bank distance from fluvial erosion. This function uses an empirical
//equation to relate specific stream power to wall or bank shear stress. It then uses an excess
//shear stress equation to calculate eroded distance.

//INPUTS:
//omega: specific stream power (W/m^2)
//width: channel bottom width (m)
//tau_c: critical shear stress of the bank material (Pa)
//bank_k: erodibility of the bank material (m^3/N/s)
//dt: time step (seconds)

//OUTPUTS:
//E: eroded distance (m)

double BankShear(const double& omega, const double& width, const double& tau_c, double& bank_k, const int& dt,
	const double& Rc, double& fluvial_factor, double& k_factor) {
	double tau_w; //wall shear stress (Pa)
	double K; //erodibility coefficient (m^3/N-s)
	double E; //Eroded distance (m)

	//empirical relationship to get wall shear from stream power - multiply by 1.5 to get max
	//value (shear at base of bank). This is based on a triangular distribution with a mean of
	//tau_w and a min of 0.
	tau_w = 0.83 * pow(omega, 0.65) * 1.5 * fluvial_factor; //Newer eq. with smooth and rough flumes

	//Adjust tau_w based on channel curvature
	if (Rc != 0) {
		if ((Rc / width >= 2) & (Rc / width < 7)) {
			//Increase calculated tau_w by factor based on Rc/w if Rc is positive (outside of bend)
			//and Rc / width < 7 (maximum value where equation predicts higher shear).
			tau_w = tau_w * 2.65 / sqrt(Rc / width);
		}
		else if ((Rc / width <= -2) & (abs(Rc) / width < 7)) {
			//Similarly, reduce wall shear stress by factor based on Rc/w if Rc is negative (inside of bend)
			tau_w = tau_w / (2.65 / sqrt(abs(Rc) / width));
		}
		else {
			tau_w = 0;
		}
	}

	//empirical relationship to get K from tau_c (from Simon et al., 2010), unless it
	//is provided by the user
	if (bank_k == 0) {
		K = 1.6 * pow(tau_c, -0.826) * 1e-6 * k_factor;
	}
	else {
		K = bank_k;
	}

	//If shear stress exceeds the critical value, calculate erosion, otherwise set erosion to zero
	if (tau_w > tau_c) {
		E = K * dt * (tau_w - tau_c);
	}
	else {
		E = 0.0;
	}

	return(E);
		
}

//BankFailure: Controls the bank failure module. Calls the Bank Stability and Toe Erosion Model (BSTEM) and
//updates bank geometry if a failure occurs.

//INPUTS:
//height: bank height (m)
//angle: bank angle (radians)
//toe_height: bank toe height (m)
//toe_angle: bank toe angle (radians)
//cohesion: bank soil cohesion (kPa)
//phi: bank friction angle (radians)
//weight: bank soil weight (kN)
//cohesion_toe: bank toe soil cohesion (kPa)
//phi_toe: bank toe friction angle (radians)
//weight_toe: bank toe soil weight (kN)
//top_width: channel top width (m)
//bottom_width: channel bottom width (m)
//fp_angle: floodplain angle (radians)
//fp_width: floodplain width (m)
//bank_tank: see "BankErosion" inputs

//OUTPUTS:
//none

void BankFailure(double& height, double& angle, double& toe_height, double& toe_angle, 
	double& cohesion, double& phi, double& weight, double& cohesion_toe, double& phi_toe,
	double& weight_toe, double& top_width, double& bottom_width, double& fp_angle, double& fp_width,
	double& bank_tank) {

	double failblock_area = 0;

	//Calculate toe_length from toe_height
	double toe_length = toe_height / sin(toe_angle);

	//Set bank geometry array (x, z)
	vector<vector<double>> bank_geom(2, vector<double>(23));
	set_bank_geom(bank_geom, height, angle, toe_length, toe_angle);

	//Set intersection points of different bank layers (just main bank and toe)
	vector<double> bank_layer_Z = { bank_geom[1][1], bank_geom[1][22], bank_geom[1][22] };
	vector<double> bank_layer_X = { bank_geom[0][1], bank_geom[0][22], bank_geom[0][22] };

	vector<double> phi_array = { phi, phi_toe };

	vector<double> cohesion_array = { cohesion, cohesion_toe };

	vector<double> unit_weight_array = { weight, weight_toe };

	//Initialize best failure plane and call BSTEM stability function
	double best_FoS = 99999, best_failure_angle = 0, best_base_Z = 0;
	calc_min_FoS(bank_geom, bank_layer_Z, bank_layer_X, cohesion_array, unit_weight_array, phi_array,
		best_FoS, best_failure_angle, best_base_Z);

	//If bank fails, update geometry and calculate eroded area
	if (best_FoS < 1.0) {
		//Determine if failure plane base is nearest bottom of bank or bottom of toe
		if (abs(best_base_Z - bank_geom[1][16]) < abs(best_base_Z - bank_geom[1][21])) {
			//Failure plane bank intersects bottom of bank (top of toe)
			//First, adjust bank x and y coordinates to incorporate fp width and angle
			bank_geom[0][0] = bank_geom[0][1] - fp_width;
			bank_geom[1][0] += fp_width * tan(fp_angle);

			//Construct polygon with coordinates
			vector<vector<double>> polygon(2, vector<double>(4));
			polygon[0][0] = bank_geom[0][16];
			polygon[1][0] = bank_geom[1][16];
			//Get coordinates of intersect of failure plane with floodplain
			int FP_point = 0;
			vector<double> FP_intersect = calc_FP_shear_intersect(bank_geom, polygon[1][0], polygon[0][0],
				best_failure_angle, FP_point);
			polygon[0][1] = FP_intersect[0];
			polygon[1][1] = FP_intersect[1];
			polygon[0][2] = bank_geom[0][1];
			polygon[1][2] = bank_geom[1][1];
			polygon[0][3] = polygon[0][0];
			polygon[1][3] = polygon[1][0];

			int num_vertices = 4;
			failblock_area = polygon_area(polygon, num_vertices);

			//Update angle
			angle = best_failure_angle;

			//Update top width
			top_width = top_width + (bank_geom[0][1] - FP_intersect[0]);

			//Update bank height and fp width
			height += (bank_geom[0][1] - FP_intersect[0]) * tan(fp_angle);
			fp_width -= (bank_geom[0][1] - FP_intersect[0]);
			
			//Don't let fp width become negative
			fp_width = max(fp_width, 0.0);
		}
		else {
			//FP intersects base of bank toe
			//First, adjust bank x and y coordinates to incorporate fp width and angle
			bank_geom[0][0] = bank_geom[0][1] - fp_width;
			bank_geom[1][0] += fp_width * tan(fp_angle);
			//Construct polygon with coordinates
			vector<vector<double>> polygon(2, vector<double>(5));
			polygon[0][0] = bank_geom[0][21];
			polygon[1][0] = bank_geom[1][21];
			//Get coordinates of intersect of failure plane with floodplain
			int FP_point = 0;
			vector<double> FP_intersect = calc_FP_shear_intersect(bank_geom, polygon[1][0], polygon[0][0],
				best_failure_angle, FP_point);
			polygon[0][1] = FP_intersect[0];
			polygon[1][1] = FP_intersect[1];
			polygon[0][2] = bank_geom[0][1];
			polygon[1][2] = bank_geom[1][1];
			polygon[0][3] = bank_geom[0][16];
			polygon[1][3] = bank_geom[1][16];
			polygon[0][4] = polygon[0][0];
			polygon[1][4] = polygon[1][0];

			int num_vertices = 5;
			failblock_area = polygon_area(polygon, num_vertices);

			//Update angle and toe angle
			angle = best_failure_angle;
			toe_angle = best_failure_angle;

			//Update top width
			top_width = top_width + (bank_geom[0][1] - FP_intersect[0]);

			//Update bank height and fp width
			height += (bank_geom[0][1] - FP_intersect[0]) * tan(fp_angle);
			fp_width -= (bank_geom[0][1] - FP_intersect[0]);

			//Don't let fp width become negative
			fp_width = max(fp_width, 0.0);
		}
	}

	//If bank failed, put failed soil block at bank toe, adjusting toe angle and bottom width
	if (failblock_area > 0) {
		FailureBlock(toe_angle, bottom_width, toe_height, failblock_area, bank_tank);
	}

}

//calc_eroded_area: Calculates the area of the bank eroded by fluvial erosion. Also calculates the area of
//any failed bank material (e.g. undercut toe causes bank collapse).

//INPUTS:
//height: bank height (m)
//toe_height: bank toe height (m)
//angle: bank angle (radians)
//toe_angle: bank toe angle (radians)
//E: eroded distance from fluvial erosion (m)
//bottom_width: bottom width of channel (m)
//top_width: top width of channel (m)
//fp_width: floodplain width (m)
//fp_angle:floodplain angle (radians)

//OUTPUTS:
//failblock_area: area of failure block (m^2)

double calc_eroded_area(double& height, double& toe_height, double& angle,
	double& toe_angle, double& E, double& bottom_width,
	double& top_width, double& fp_width, double& fp_angle) {

	//Calculate different components of geometry of eroded area
	double y = toe_height / tan(toe_angle);
	double x = E - y;
	double z = (height - toe_height) / sin(angle);
	double w = sqrt(x * x + z * z - 2 * x * z * cos(angle));

	//Calculate new bank angle
	double new_angle = angle + asin(x * sin(angle) / w);
		 
	//Calculate failed bank eroded area
	double failblock_area = 0.5 * x * toe_height + 0.5 * (height - toe_height) * x;

	//Check if new_angle is greater than 90 deg. If so, update area of failure block and new_angle
	if (new_angle > M_PI / 2.0) {
		double x2 = x - sqrt(z * z - ((height - toe_height) * (height - toe_height)));
		failblock_area += 0.5 * x2 * (height - toe_height + x2 * tan(fp_angle));
		new_angle = M_PI / 2.0;
		//update top width, fp width, and bank height
		top_width += x2;
		height += x2 * tan(fp_angle);
		fp_width -= x2;
		//Dont let fp_width become negative
		fp_width = max(fp_width, 0.0);
	}

	//Set bank angle to new_angle
	angle = new_angle;

	return(failblock_area);
}

//FluvialErosion: Calculates fluvial erosion and updates channel geometry accordingly.

//INPUTS:
//omega: specific stream power (W/m^2)
//toe_angle: angle of bank toe (radians)
//angle: bank angle (radians)
//toe_height: height of bank toe (m)
//height: bank height (m)
//eroded_area: area of bank sediment eroded by fluvial erosion (m^2)
//bottom_width: channel bottom width (m)
//tau_c: critical shear stress of bank toe material (Pa)
//bank_k: erodibility of bank toe material (m^3/N/s)
//dt: time step (s)
//Rc: radius of curvature of channel bends (m)
//top_width: channel top width (m)
//fp_width: floodplain width (m)
//fp_angle: floodplain angle (radians)
//bank_tank: see inputs to "BankErosion"
//skip_bank_erosion: see inputs to "BankErosion"

//OUTPUTS:
//E: bank retreat from fluvial erosion (m)

double FluvialErosion(double& omega, double& toe_angle, double& angle, double& toe_height, double& height,
	double& eroded_area, double& bottom_width, double& tau_c, double& bank_k, int& dt, double& Rc,
	double& top_width, double& fp_width, double& fp_angle, double& bank_tank, int& skip_bank_erosion,
	double& fluvial_factor, double& k_factor) {

	//Get eroded distance from bank shear function
	double E = BankShear(omega, bottom_width, tau_c, bank_k, dt, Rc, fluvial_factor, k_factor);

	//Erode material from "bank_tank", if there is any
	double eroded_area_tank = 0;
	if (bank_tank > 0) {
		//Assume rectangle with a height of "toe_height"
		eroded_area_tank = E * toe_height;

		//Don't let eroded_area_tank exceed tank capacity
		if (eroded_area_tank > bank_tank) {
			//Update eroded distance
			E -= bank_tank / toe_height;
			eroded_area_tank = bank_tank;
			bank_tank = 0;
		}
		else {
			bank_tank -= eroded_area_tank;
		}
	} 

	//Don't let E be larger than the distance from the bank toe to the edge of fp (if fp is zero)
	if (fp_width == 0) {
		E = min(E, (height - toe_height) / tan(angle) + toe_height / tan(toe_angle));
		skip_bank_erosion = 2;
	}

	//Only update geometry if bank_tank is zero
	if (bank_tank == 0) {
		//Update bottom width based on erosion of single bank
		bottom_width = bottom_width + E;

		//Update toe angle based on erosion magnitude - store old toe_angle
		double old_toe_angle = toe_angle;
		toe_angle = atan(toe_height / (toe_height / tan(toe_angle) - E));

		//Check if toe angle is greater than 90 degrees. If it is, erode bank toe point to give a toe angle of
		//90, and update bank angle
		if ((toe_angle > M_PI / 2.0) | (toe_angle < 0)) {
			//Update bank angle and get eroded area
			double failblock_area = calc_eroded_area(height, toe_height, angle,
				old_toe_angle, E, bottom_width, top_width, fp_width, fp_angle);

			//Update toe bank angle
			toe_angle = M_PI / 2.0;

			//Deposit failed bank material at toe and update toe angle and width
			FailureBlock(toe_angle, bottom_width, toe_height, failblock_area, bank_tank);
		}

		//Calculate eroded area of toe (triangle)
		eroded_area = 0.5 * toe_height * E;
	}

	//Add any eroded material from the bank tank to the total eroded area
	eroded_area += eroded_area_tank;

	return(E);
}

//FailureBlock: Calculates the size of the soil block from bank failure (or bank collapse from an undercut toe) and
//updates toe geometry and channel width accordingly.

//INPUTS:
//toe_angle: angle of bank toe (radians)
//bottom_width: channel bottom width (m)
//toe_height: height of bank toe (m)
//failblock_area: area of failure block (m^2)
//bank_tank: see "BankErosion" inputs

//OUTPUTS: none

void FailureBlock(double& toe_angle, double& bottom_width, double& toe_height, double& failblock_area, double& bank_tank) {

	double eroded_dist = 2.0 * failblock_area / toe_height;

	//Don't let eroded distance exceed 95% of the channel half-width or (channel half-width - 10 cm), whichever
	//is smaller
	if (eroded_dist > min(bottom_width / 2.0 - 0.1, 0.95 * bottom_width / 2.0)) {
		eroded_dist = min(bottom_width / 2.0 - 0.1, 0.95 * bottom_width / 2.0);

		//Add extra sediment to bank_tank
		bank_tank += failblock_area - 0.5 * eroded_dist * toe_height;
	}

	//Update toe angle and bottom width
	toe_angle = atan(toe_height / (toe_height / tan(toe_angle) +
		eroded_dist));
	bottom_width -= eroded_dist; //Shrink bottom width

}

//BankErosion calculates bank erosion at each cross section from both fluvial erosion and mass
//wasting. Channel width is updated based on simulated erosion.

//INPUTS:
//tau_c: critical shear stress of the bank material (Pa) (by reach)
//bank_k: erodibility of the bank material (m^3/N/s) (by reach)
//q: vector of unit discharge (m^2/s) at each cross section
//bottom_width: vector of channel bottom width by each simulation day and cross section (m)
//top_width: vector of channel top width by each simulation day and cross section (m)
//slope: vector of channel slope at each cross section (m/m)
//dt: time step (s)
//n_nodes: number of channel nodes/reaches
//i: model iteration
//height_RB: right bank heights by XS (m)
//height_LB: left bank heights by XS (m)
//angle_RB: right bank angles by XS (radians)
//angle_LB: left bank angles by XS (radians)
//toe_height_RB: right bank toe heights by XS (m)
//toe_height_LB: left bank toe heights by XS (m)
//toe_angle_RB: right bank toe angles by XS (radians)
//toe_angle_LB: left bank toe angles by XS (radians)
//cohesion: bank cohesion by reach (kPa)
//phi: bank friction angle by reach (radians)
//weight: bank soil weight by reach (kN)
//cohesion_toe: bank toe cohesion by reach (kPa)
//phi_toe: bank toe friction angle by reach (radians)
//weight_toe: bank toe soil weight by reach (kN)
//bank_erosion: type of bank erosion to model ("none", "fluvial", "failure", or "both")
//n_xs: number of cross sections to be modeled (for each reach)
//bank_sed: the amount of sediment from bank erosion to be added to bed material load (m^3/s)
//dx_array: XS spacing (m)
//Rc: channel bend radius of curvature by reach (m)
//sinuosity: channel sinuosity by reach
//n_bends: number of channel bends per reach
//meandering: whether to include meandering in the model (1) or not (0)
//bank_bedload_prop: proportion of the bank sediment to be added to bedload, by reach (0-1)
//day: day of the simulation
//day_old: day of previous time step
//bank_sed_mass: mass of washload from bank erosion (kg)
//bank_p_mass: mass of phosphorus (or other pollutant) loading from bank erosion (kg)
//p_conc: concentration of phosphorus (or other pollutant) in bank material (mg/kg soil)
//LB_x: station of the bottom of the left bank toe (for plotting XS)
//fp_angle: angle of floodplain by reach (radians)
//fp_width_R: width of right floodplain by XS (m)
//fp_width_L: width of left floodplain by XS (m)
//bank_tank_RB: storage "tank" of failed bank material for right bank. Used to conserve mass (m^2)
//bank_tank_LB: same as above but for left bank
//skip_bank_erosion_L/R: whether to skip modeling fluvial (1) or both types (2) of bank erosion. Used to 
		//not erode channel beyond edge of floodplain.
//lambda: porosity of bed material (assumed to be the same as bedload material in bank)
//bank_sed_vol: volume of eroded washload (m^3)

//OUTPUTS: none

void BankErosion(vector<double>& tau_c, vector<double>& bank_k, vector<vector<double>>& q, 
	vector<vector<double>>& bottom_width, vector<vector<double>>& top_width, 
	vector<vector<double>>& slope, int& dt, const int & n_nodes, 
	const long& i,vector<vector<double>>& height_RB, vector<vector<double>>& height_LB, 
	vector<vector<double>>& angle_RB, vector<vector<double>>& angle_LB, vector<vector<double>>&toe_height_RB, 
	vector<vector<double>>&toe_height_LB, vector<vector<double>>& toe_angle_RB, 
	vector<vector<double>>& toe_angle_LB, vector<double>& cohesion, vector<double>& phi, 
	vector<double>& weight, vector<double>& cohesion_toe, vector<double>& phi_toe, 
	vector<double>& weight_toe,	string& bank_erosion, vector<int>& n_xs, 
	vector<vector<double>>& bank_sed, vector<vector<double>>& dx_array, 
	vector<vector<double>>& Rc, vector<double>& sinuosity, 
	vector<vector<double>>& n_bends, int& meandering, vector<double>& bank_bedload_prop,
	int& day, int& day_old, vector<vector<double>>& bank_sed_mass, vector<vector<double>>& bank_p_mass,
	vector<double>& p_conc, vector<vector<double>>& LB_x, vector<vector<double>>& fp_angle,
	vector<vector<double>>& fp_width_R, vector<vector<double>>& fp_width_L, vector<vector<double>>& bank_tank_RB,
	vector<vector<double>>& bank_tank_LB, vector<vector<int>>& skip_bank_erosion_L, 
	vector<vector<int>>& skip_bank_erosion_R, double& lambda, double& bank_sed_vol, double& fluvial_factor,
	double& k_factor, vector<vector<double>>& meander_erosion) {

	//Initialize stream power, eroded distance, and eroded area to zero
	double omega = 0;
	double E = 0, eroded_area = 0;
	double Rc_RB = 0, Rc_LB = 0; //Radius of curvature of right and left banks

	//For each cross section, calculate stream power and call BankShear function
	for (int j = 0; j < n_nodes; ++j) {
		for (int k = 0; k < n_xs[j]; ++k) {
			//Initialize eroded_area to zero
			eroded_area = 0;
			double eroded_area_LB = 0, eroded_area_RB = 0;
			double eroded_area_meandering = 0;

			//Fluvial erosion
			if (bank_erosion != "failure") {
				omega = slope[j][k] * q[j][k] * 9810.0;
				Rc_RB = 0;
				Rc_LB = 0;

				//Right bank erosion
				if (skip_bank_erosion_R[j][k] != 2) {
					double E_RB = FluvialErosion(omega, toe_angle_RB[j][k], angle_RB[j][k], toe_height_RB[j][k],
						height_RB[j][k], eroded_area_RB, bottom_width[j][k], tau_c[j], bank_k[j], dt, Rc_RB,
						top_width[j][k], fp_width_R[j][k], fp_angle[j][k], bank_tank_RB[j][k], skip_bank_erosion_R[j][k],
						fluvial_factor, k_factor);
				}

				//Left bank erosion
				if (skip_bank_erosion_L[j][k] != 2) {
					double width_old = bottom_width[j][k];
					double E_LB = FluvialErosion(omega, toe_angle_LB[j][k], angle_LB[j][k], toe_height_LB[j][k],
						height_LB[j][k], eroded_area_LB, bottom_width[j][k], tau_c[j], bank_k[j], dt, Rc_LB,
						top_width[j][k], fp_width_L[j][k], fp_angle[j][k], bank_tank_LB[j][k], skip_bank_erosion_L[j][k],
						fluvial_factor, k_factor);
					LB_x[j][k] = LB_x[j][k] - (bottom_width[j][k] - width_old);
				}
				
				if (toe_angle_LB[j][k] != toe_angle_LB[j][k] | (toe_angle_RB[j][k] != toe_angle_RB[j][k])) {
					cout << "Toe angle banks\n";
				}
				//Update Rc, sinuosity, and dx based on bend erosion
				if ((meandering == 1) & (sinuosity[j] <= 2.5) & 
					(skip_bank_erosion_L[j][k] != 2) & (skip_bank_erosion_R[j][k] != 2)) {
					//Outer bank erosion distance
					double E_outer = BankShear(omega, bottom_width[j][k], tau_c[j], bank_k[j],
						dt, Rc[j][k], fluvial_factor, k_factor);
					meander_erosion[j][k] += E_outer;
					////Inner bank erosion distance
					//double E_inner = BankShear(omega, bottom_width[j][k], tau_c[j], bank_k[j],
					//	dt, -Rc[j][k]);

					//double E_avg = (E_outer - E_inner) / 2.0;
					double E_avg = E_outer;
					if ((E_avg > 0)){ 

						bool angle_adj = false;
						//Calculate angle of circular bend
						double bend_angle = dx_array[j][k] / n_bends[j][k] / Rc[j][k];
						//Calculate segment length
						double seg_length = 2.0 * Rc[j][k] * sin(0.5 * bend_angle);
						//adjust angle if angle greater than 180 deg
						if (bend_angle > M_PI) { angle_adj = true; }

						//Calculate initial area of curve
						double A_initial = 0.5 * Rc[j][k] * Rc[j][k] * (bend_angle - sin(bend_angle));

						//Calculate bend height
						double h = Rc[j][k] - 0.5 * sqrt(4 * Rc[j][k] * Rc[j][k] - seg_length * seg_length);
						//Adjust h depending on the angle
						if (angle_adj) {
							h -= E_avg;
						}
						else {
							h += E_avg;
						}
						//Calculate new radius of curvature
						double Rc_old = Rc[j][k];
						Rc[j][k] = (4 * h * h + seg_length * seg_length) / (8 * h);

						//Calculate new bend angle
						double angle_old = bend_angle;
						bend_angle = 2 * asin(seg_length / (2 * Rc[j][k])); 

						//Adjust angle if original was greater than 180 deg, or if new angle is less than old angle
						if (angle_adj | (bend_angle < angle_old)) { bend_angle = 2.0 * M_PI - bend_angle; }

						if (bend_angle < angle_old) {
							cout << "MEANDERING!\n";
						}

						//Calculate final area of curve
						double A_final = 0.5 * Rc[j][k] * Rc[j][k] * (bend_angle - sin(bend_angle));

						//Calculate eroded area from meandering (outer bend erosion - inner bend deposition)
						double A_outer = meandering_eroded_area(dx_array[j][k], n_bends[j][k],
							Rc[j][k] + bottom_width[j][k] / 2.0, E_avg);
						double A_inner = meandering_eroded_area(dx_array[j][k], n_bends[j][k],
							Rc[j][k] - bottom_width[j][k] / 2.0, E_avg);

						//eroded_area_meandering = (A_outer - A_inner) * n_bends[j][k] * 
						//	(height_RB[j][k] + height_LB[j][k]) / 2.0;

						//eroded_area_meandering += (A_final - A_initial) * n_bends[j][k] / 2.0 * height_RB[j][k];
						//eroded_area_meandering += (A_final - A_initial) * n_bends[j][k] / 2.0 * height_LB[j][k];

						if (eroded_area_meandering < 0) {
							eroded_area_meandering = 0;
							//cout << "MEANDERING!\n";
						}
						//Calculate arc length
						double arc_length = Rc[j][k] * bend_angle;

						//Update dx with the new arc length times the number of bends
						double dx_old = dx_array[j][k];
						dx_array[j][k] = arc_length * n_bends[j][k];

						//area
						vector<double> chnl_x(10);
						vector<double> chnl_y(10);
						double z = 0.0;
						double depth = (height_LB[j][k] + height_RB[j][k]) / 2.0;
						double area_fp, area_chnl, perim, perim_fp, perim_chnl, ws_width;
						create_chnl_geom(height_LB[j][k], toe_height_LB[j][k], angle_LB[j][k], toe_angle_LB[j][k],
							height_RB[j][k], toe_height_RB[j][k], angle_RB[j][k], toe_angle_RB[j][k], bottom_width[j][k],
							fp_angle[j][k], fp_width_R[j][k], fp_width_L[j][k], z, chnl_x, chnl_y);
						double chnl_area = calc_flow_area(chnl_x, chnl_y, depth,
							area_fp, area_chnl, perim, perim_fp, perim_chnl, ws_width);

						eroded_area_meandering = area_chnl * (dx_array[j][k] - dx_old);
					}
				}
			}

			//Check bank stability
			if ((bank_erosion != "fluvial") & (day != day_old)) {
				if (skip_bank_erosion_R[j][k] == 0 & fp_width_R[j][k] != 0) {
					//run BSTEM for right bank
					BankFailure(height_RB[j][k], angle_RB[j][k], toe_height_RB[j][k],
						toe_angle_RB[j][k], cohesion[j], phi[j], weight[j], cohesion_toe[j],
						phi_toe[j], weight_toe[j], top_width[j][k], bottom_width[j][k], fp_angle[j][k],
						fp_width_R[j][k], bank_tank_RB[j][k]);
				}

				//run BSTEM for left bank
				if (skip_bank_erosion_L[j][k] == 0 & fp_width_L[j][k] != 0) {
					double width_old = bottom_width[j][k];
					BankFailure(height_LB[j][k], angle_LB[j][k], toe_height_LB[j][k],
						toe_angle_LB[j][k], cohesion[j], phi[j], weight[j], cohesion_toe[j],
						phi_toe[j], weight_toe[j], top_width[j][k], bottom_width[j][k], fp_angle[j][k],
						fp_width_L[j][k], bank_tank_LB[j][k]);
					LB_x[j][k] = LB_x[j][k] - (bottom_width[j][k] - width_old);
				}
			}

			//Calculate bed material loading rate from eroded bank material, converting area to
			//volume, dividing by s/day, giving final units of m^3/s.
			bank_sed[j][k] = ((eroded_area_RB + eroded_area_LB) * dx_array[j][k] + eroded_area_meandering) * bank_bedload_prop[j]  * (1 - lambda) / dt;

			//Calculate sediment and phosphorus mass loading from bank erosion
			bank_sed_vol += ((eroded_area_RB + eroded_area_LB) * dx_array[j][k] + eroded_area_meandering) * (1 - bank_bedload_prop[j]); //m^3 sediment
			bank_sed_mass[j][k] = ((eroded_area_RB + eroded_area_LB) * dx_array[j][k] + eroded_area_meandering) * (1 - bank_bedload_prop[j]) *
				1600; //kg sediment; assumes soil bulk density of 1.6 Mg/kg
			bank_p_mass[j][k] = bank_sed_mass[j][k] * p_conc[j] / 1000000.0; //kg phosphorus
		}

	}

	//Loop through and update sinuosity
	if (meandering == 1) {
		for (int j = 0; j < n_nodes; j++) {
			double Lv = 0; //Valley length
			double Ls = 0; //Stream length
			for (int k = 0; k < n_xs[j]; k++) {
				double bend_angle = dx_array[j][k] / n_bends[j][k] / Rc[j][k];
				Lv += 2 * Rc[j][k] * sin(0.5 * bend_angle) * n_bends[j][k];
				Ls += dx_array[j][k];
			}
			double sinuosity_old = sinuosity[j];
			sinuosity[j] = Ls / Lv;
			//if (sinuosity[j] < sinuosity_old) {
			//	cout << "SINUOSITY!\n";
			//}
		}
	}

}

