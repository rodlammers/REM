#include "stdafx.h"
using namespace std;

void set_bank_geom(vector<vector<double>>& bank_geom, double& height, double& angle,
	double& toe_length, double& toe_angle) {

	//Set initial bank geometry based on height and angle (1st row x, 2nd row z)
	bank_geom[0][0] = 0;
	bank_geom[0][1] = height;
	bank_geom[1][0] = height;
	bank_geom[1][1] = height - 0.000001;
	bank_geom[1][21] = 0;
	bank_geom[1][22] = -0.000001;

	//Set toe z coordinates
	for (int i = 0; i < 5; ++i) {
		bank_geom[1][16 + i] = (5 - i) / 5.0 * sin(toe_angle) * toe_length;
	}

	//Set main bank coordinates
	for (int i = 0; i < 15; ++i) {
		bank_geom[0][i + 2] = bank_geom[0][1] + (i + 1) / 15.0 * (bank_geom[1][1] - bank_geom[1][16]) / tan(angle);
		if (i != 14) {
			bank_geom[1][i + 2] = bank_geom[1][16] + ((14 - i) / 15.0 * (bank_geom[1][1] - bank_geom[1][16]));
		}
	}

	//Set toe x coordinates
	for (int i = 0; i < 5; ++i) {
		bank_geom[0][17 + i] = bank_geom[0][16] + (i + 1) / 5.0 * cos(toe_angle) * toe_length;
	}
	bank_geom[0][22] = bank_geom[0][21] + 1;
}

double calc_min_angle(double& base_Z, vector<double>& Phi, vector<vector<double>>& bank_geom,
	vector<double>& bank_layer_Z, vector<double>& bank_layer_X) {

	double min_angle = 0, length = 0, sum_length = 0;

	double num_layers = calc_num_layers(bank_geom, base_Z, bank_layer_Z);

	for (int i = 0; i < num_layers - 1; ++i) {
		length = bank_layer_Z[i] - bank_layer_Z[i + 1];
		sum_length += length;
		min_angle += Phi[i] * length;
	}
	length = bank_layer_Z[num_layers - 1] - base_Z;
	sum_length += length;
	min_angle += Phi[num_layers - 1] * length;
	if (sum_length == 0) {
		min_angle = Phi[num_layers - 1];
	}
	else {
		min_angle = min_angle / sum_length;
	}

	//Adjust by 0.75 because failure surface could be smaller than the friction angle for saturated conditions
	min_angle = 0.75 * min_angle;

	if (M_PI / 90.0 > min_angle) { min_angle = M_PI / 90; }

	return(min_angle);
}

double calc_max_angle(double& base_Z, double& base_X, vector<vector<double>>& bank_geom) {

	//Find first node above the failure plane base
	int point = 22;
	int limit_point = 0;
	double max_angle = 0;

	while (point > 0) {
		if ((base_Z >= bank_geom[1][point]) & (bank_geom[1][point - 1] > base_Z)) {
			limit_point = point - 1;
		}
		point -= 1;
	}

	double length = base_X - bank_geom[0][limit_point];

	if (0 >= length) {
		max_angle = 0.5 * M_PI;
	}
	else {
		max_angle = atan((bank_geom[1][limit_point] - base_Z) / length);
		for (int i = 1; i < limit_point - 1; ++i) {
			length = base_X - bank_geom[0][i];
			if (length > 0) {
				max_angle = min(max_angle, atan((bank_geom[1][i] - base_Z) / length));
			}
		}
	}

	if (0 >= max_angle) { max_angle = 0.5 * M_PI; }

	//Reduce max angle: minimum angle difference is 4 degrees
	double reduced_angle = 2 * M_PI / 90;
	if (0.5 * M_PI > max_angle) { max_angle = max_angle - reduced_angle; }

	return(max_angle);
}

double set_bank_intersect_X(vector<vector<double>>& bank_geom, double& base_Z) {
	double base_X = 0;
	int point = 22;
	double Z1 = 0, Z2 = 0, X1 = 0, X2 = 0;

	while (point > 0) {
		if ((base_Z > bank_geom[1][point]) & (bank_geom[1][point - 1] >= base_Z)) {
			Z1 = bank_geom[1][point - 1];
			Z2 = bank_geom[1][point];
			X1 = bank_geom[0][point - 1];
			X2 = bank_geom[0][point];
		}
		point = point - 1;
	}

	base_X = X1 + (X2 - X1) * (base_Z - Z1) / (Z2 - Z1);

	return(base_X);
}

vector<double> calc_pore_water_force(double& base_Z, double& failure_angle,
	vector<double>& Phi, double& Phib, vector<vector<double>>& bank_geom, double& FP_fail_Z,
	vector<double>& bank_layer_Z, vector<double>& unit_weight) {

	double mid_Z = 0, length = 0;
	double angular_coeff = -0.063;
	vector<double> pore_water_force(2);
	double GW_Z = bank_geom[1][1] / 2.0; //GW elevation is half of total bank height

										 //For main bank layer and toe layer
	for (int i = 0; i < 2; ++i) {
		//If failure block includes layer, set pore_force, else set to zero
		if (bank_layer_Z[i] > base_Z) {
			mid_Z = bank_layer_Z[i] - 0.5 * (bank_layer_Z[i] - max(base_Z, bank_layer_Z[i + 1]));
			pore_water_force[i] = -9.807 * (GW_Z - mid_Z);
		}
		else {
			pore_water_force[i] = 0;
		}

		//Adjust bulk unit weight?? - code in BSTEM
		if (pore_water_force[i] > 0) { unit_weight[i] += angular_coeff * pore_water_force[i]; }
		if (0 > unit_weight[i]) { unit_weight[i] = 0; }

		//Compute length of failure plane within this layer
		if (FP_fail_Z >= bank_layer_Z[i]) {
			length = (bank_layer_Z[i] - max(base_Z, bank_layer_Z[i + 1])) / sin(failure_angle);
		}
		else {
			length = (FP_fail_Z - bank_layer_Z[i + 1]) / sin(failure_angle);
		}
		if (0 > length) { length = 0; }

		//Comput matric suction or positive PWP
		if (pore_water_force[i] > 0) {
			pore_water_force[i] = -length * pore_water_force[i] * tan(Phib);
		}
		else {
			pore_water_force[i] = -length * pore_water_force[i] * tan(Phi[i]);
		}
	}

	return(pore_water_force);
}

double set_water_bank_intersect(vector<vector<double>>& bank_geom, double& WS_Z) {

	//Go to row with top of the bank and search for points in bank profile surrounding
	//the water surface elevation
	double Z1 = 0, Z2 = 0, X1 = 0, X2 = 0, WS_X = 0;

	int point = 0;
	while (bank_geom[1].size() > (point + 1)) {
		Z1 = bank_geom[1][point];
		Z2 = bank_geom[1][point + 1];
		if (WS_Z >= Z2) { break; }
		point += 1;
	}
	X1 = bank_geom[0][point];
	X2 = bank_geom[0][point + 1];

	//Check for bankfull or over-bankfull
	if (WS_Z > Z1) {
		WS_X = bank_geom[0][1];
	}

	//Calculate station of intersect point
	WS_X = X1 + (X2 - X1) * (WS_Z - Z1) / (Z2 - Z1);

	return(WS_X);
}

vector<double> calc_FP_shear_intersect(vector<vector<double>>& bank_geom, double& base_Z, double& base_X,
	double& failure_angle, int& FP_point) {

	double distance = 0, Z1 = 0, Z2 = 0, X1 = 0, X2 = 0, tan_angle = 0, tan_failure_angle = 0, tan_FP_angle = 0, top_Z = 0, top_X = 0;
	FP_point = 21;
	//Set tangent of failure surface
	tan_failure_angle = tan(failure_angle);

	//Do calculation depending on whether failure angle is 90 degrees
	if (0.5 * M_PI > failure_angle) {


		//Loop to find first point above failure plane base point
		while (FP_point > 0) {
			if ((bank_geom[1][FP_point] >= base_Z) & (base_X > bank_geom[0][FP_point])) { break; }
			FP_point -= 1;
		}

		//Now loop through remaining points to find where top of failure plane intersect bank or fp
		while (FP_point > 0) {
			distance = base_X - bank_geom[0][FP_point];
			if (abs(distance) < 0.00001) { distance = 0.00001; }
			tan_angle = (bank_geom[1][FP_point] - base_Z) / distance;
			if (0 >= tan_angle) { tan_angle = 999999; }
			if (tan_failure_angle > tan_angle) { break; }
			FP_point -= 1;
		}

		//Set elevation and stations for points bounding FP_point
		X1 = bank_geom[0][FP_point];
		X2 = bank_geom[0][FP_point + 1];
		Z1 = bank_geom[1][FP_point];
		Z2 = bank_geom[1][FP_point + 1];

		tan_FP_angle = (Z1 - Z2) / max(0.00001, X2 - X1);

		if ((tan_FP_angle >= tan_failure_angle) | (0 > tan_FP_angle)) {
			top_Z = bank_geom[1][FP_point];
			top_X = base_X - (top_Z - base_Z) / tan_failure_angle;
		}
		else {
			top_X = (base_Z - Z2 + base_X * tan_failure_angle - X2 * tan_FP_angle) /
				(tan_failure_angle - tan_FP_angle);
			top_Z = base_Z - tan_failure_angle * (top_X - base_X);
		}

		if (bank_geom[0][0] > top_X) {
			bank_geom[0][0] = top_X;
			bank_geom[1][0] = top_Z;
		}
	}
	else {
		//Loop up bank points to find first node above the failure plane base
		int FP_point = 21;
		while (FP_point > 1) {
			if ((bank_geom[1][FP_point] > base_Z) & (base_Z >= bank_geom[1][FP_point + 1])) { break; }
			FP_point -= 1;
		}

		//Continue looping up the bank until the first node with a station less than the base station is found
		while (FP_point > 1) {
			X1 = bank_geom[0][FP_point];
			X2 = bank_geom[0][FP_point - 1];
			if (base_X > X2) { break; }
			FP_point -= 1;
		}
		Z1 = bank_geom[1][FP_point];
		Z2 = bank_geom[1][FP_point - 1];

		//Calculation station and elevation of intersection point
		top_X = base_X;
		top_Z = Z1 + (Z2 - Z1) * (top_X - X1) / (X2 - X1);
	}

	//Combine top station and elevation into a vector and return
	vector<double> FP_intersect = { top_X, top_Z };

	return(FP_intersect);
}

vector<double> calc_weight(vector<vector<double>>& bank_geom, vector<double>& FP_intersect, double& base_X,
	double& base_Z, vector<double>& unit_weight, vector<double>& bank_layer_Z, vector<double>& bank_layer_X,
	double& failure_angle, int& FP_point) {

	vector<double> weight = { 0, 0 };
	double area = 0;
	vector<vector<double>> layer_geom(2, vector<double>(2 * bank_geom[0].size() + 4)); //Array to store bank layer geometry to calc area
	int point = 0;

	//Find first layer within the failure block
	int num_layers = calc_num_layers(bank_geom, base_Z, bank_layer_Z);
	int first_layer = 0;
	for (int layer = 0; layer < num_layers; ++layer) {
		if (bank_layer_Z[layer + 1] >= FP_intersect[1]) {
			first_layer = layer + 1;
			break;
		}
	}

	//Initialize two vertices
	double X_BL = FP_intersect[0];
	double Z_BL = FP_intersect[1];

	for (int layer = first_layer; layer < num_layers; ++layer) {
		//Top valley-side intersect of soil layer and slip surface is the bottom valley-side
		//intersect of the soil layer above
		double X_TL = X_BL;
		double Z_TL = Z_BL;

		//Assemble vertices of polygon of soil layer within the failure block

		//Bottom side
		int vertex = -1;
		if (layer == num_layers - 1) {
			//One vertex only: intersect of slip surface and bank face
			layer_geom[0][vertex + 1] = base_X;
			layer_geom[1][vertex + 1] = base_Z;
			vertex += 1;
		}
		else {
			//Two vertices, stream side:
			layer_geom[0][vertex + 1] = bank_layer_X[layer + 1];
			layer_geom[1][vertex + 1] = bank_layer_Z[layer + 1];
			//Valley side: intersect of bottom of soil layer and slip surface
			X_BL = base_X - (bank_layer_Z[layer + 1] - base_Z) / tan(failure_angle);
			Z_BL = bank_layer_Z[layer + 1];
			layer_geom[0][vertex + 2] = X_BL;
			layer_geom[1][vertex + 2] = Z_BL;
			vertex += 2;
		}

		//Top valley-side vertex
		layer_geom[0][vertex + 1] = X_TL;
		layer_geom[1][vertex + 1] = Z_TL;
		vertex += 1;

		//Top stream-side vertex
		if (layer > first_layer) {
			layer_geom[0][vertex + 1] = bank_layer_X[layer];
			layer_geom[1][vertex + 1] = bank_layer_Z[layer];
			vertex += 1;

			//Find first node on bank face that may possibly be added to the polygon
			point = 22;
			while (layer_geom[1][vertex] > bank_geom[1][point]) {
				point -= point;
			}
		}
		else {
			point = FP_point; // Point where failure plane intersects floodplain?
		}

		//Add bank nodes to polygon until a node's elevation equals or dops below the elevation 
		//of the intersect of bank profile with the bottom of soil layer.
		while (bank_geom[1][point] > layer_geom[1][0]) {

			if (layer_geom[1][vertex] < bank_geom[1][point]) {
				point += 1;
			}
			else {
				layer_geom[0][vertex + 1] = bank_geom[0][point];
				layer_geom[1][vertex + 1] = bank_geom[1][point];
				vertex += 1;
				point += 1;
			}
		}
		layer_geom[0][vertex + 1] = layer_geom[0][0];
		layer_geom[1][vertex + 1] = layer_geom[1][0];
		vertex += 1;

		//Comput area of failure block and convert to weight
		weight[layer] = polygon_area(layer_geom, vertex) * unit_weight[layer];
	}



	//Loop through main bank layer and toe and calculate weight of soil in failure block
	//if failure block includes that layer
	//for (int i = 0; i < 2; ++i) {
	//	//Is failure plane top above bottom of bank layer?
	//	if (FP_intersect[1] > bank_layer_Z[i + 1]) {
	//		//Is failure plane bottom below bottom of bank layer?
	//		if (base_Z < bank_layer_Z[i + 1]) {
	//			//Cacluate station where layer line intersects failure plane
	//			double failure_layer_int_X = FP_intersect[0] + (bank_layer_Z[i + 1] - FP_intersect[1]) *
	//				(base_X - FP_intersect[0]) / (base_Z - FP_intersect[1]);
	//			if (FP_intersect[1] >= bank_layer_Z[i]) {
	//				//If polygon then calculate area
	//				vector<vector<double>> polygon = { { FP_intersect[0], bank_geom[0][1], bank_geom[0][16],
	//				failure_layer_int_X }, { FP_intersect[1], bank_geom[1][1], bank_geom[1][16],
	//				bank_geom[1][16]} };
	//				area = polygon_area(polygon);
	//			}
	//			else {
	//				//Area of triangle
	//				area = 0.5 * (bank_geom[0][16] - failure_layer_int_X) * 
	//					(FP_intersect[1] - bank_geom[1][16]);
	//			}
	//		}
	//		else {
	//			//area of triangle
	//			area = 0.5 * (bank_geom[0][1] - FP_intersect[0]) * (FP_intersect[1] - base_Z);
	//		}

	//	}

	//	//Calculate weight
	//	weight[i] = area * unit_weight[i];
	//}

	return(weight);
}

double polygon_area(vector<vector<double>>& polygon, int& num_vertices) {

	double polygon_area = 0;

	for (int i = 0; i < num_vertices; i++) {
		int i_plus = i + 1;
		if (i_plus >= num_vertices) { i_plus = 0; }
		polygon_area = polygon_area + polygon[0][i] * polygon[1][i_plus] -
			polygon[1][i] * polygon[0][i_plus];
	}

	polygon_area = 0.5 * polygon_area;
	if (polygon_area < 0) { polygon_area = -polygon_area; }

	return(polygon_area);
}

void calc_water_force(vector<vector<double>>& bank_geom, double& base_X, double& base_Z,
	vector<double>& confining_force, vector<double>& confining_angle,
	vector<double>& bank_layer_Z, vector<double>& bank_layer_X) {

	//Assume water surface elevation is at the top of the toe - for now...
	double WS_Z = bank_geom[1][1] / 2.0;
	double WS_X = set_water_bank_intersect(bank_geom, WS_Z);

	//Check if water surface is below the base of failure plane surface and set all confining
	//forces and angles to zero.
	if (base_Z >= WS_Z) {
		for (int i = 0; i < 2; ++i) {
			confining_force[i] = 0;
			confining_angle[i] = 0;
		}
	}
	else {
		//Loop over each layer to find confining force and angle

		//First find point on bank just above failure plane base
		int point = 0;
		while (bank_geom[1].size() > point) {
			if (bank_geom[1][point] < base_Z) { break; }
			point += 1;
		}
		point -= 1;

		int num_layers = calc_num_layers(bank_geom, base_Z, bank_layer_Z); //Number of layers in failure block

																		   //Loop through layers
		for (int i = (num_layers - 1); i >= 0; i--) {
			//Initialize horizontal and veritcal force components
			double force_H = 0, force_V = 0;
			double bottom_Z = 0, bottom_X = 0, top_Z = 0, top_X = 0;

			//Set top and bottom coordinates of soil layer
			if (i == num_layers - 1) {
				//Base of failure surface
				bottom_Z = base_Z;
				bottom_X = base_X;
			}
			else {
				//Bottom of soil layer
				bottom_Z = bank_layer_Z[i + 1];
				bottom_X = bank_layer_X[i + 1];
			}

			if (WS_Z < bank_layer_Z[i]) {
				//Intersect of water surface elevation and bank profile
				top_Z = WS_Z;
				top_X = WS_X;
			}
			else {
				//Top of soil layer
				top_Z = bank_layer_Z[i];
				top_X = bank_layer_X[i];
			}

			//Compute confining force for each point on bank profile
			double Z1 = 0, Z2 = 0, X1 = 0, X2 = 0, top_pressure = 0, bottom_pressure = 0,
				delta_X = 0, delta_Z = 0, pressure = 0;
			Z1 = bottom_Z;
			X1 = bottom_X;
			top_pressure = 9.807 * (WS_Z - Z1);
			if (point > 0) {
				while (bank_geom[1][point] < top_Z) {
					Z2 = Z1;
					X2 = X1;
					bottom_pressure = top_pressure;
					Z1 = bank_geom[1][point];
					X1 = bank_geom[0][point];
					top_pressure = 9.807 * (WS_Z - Z1);
					delta_X = X2 - X1;
					delta_Z = Z1 - Z2;
					pressure = 0.5 * (bottom_pressure + top_pressure);
					force_H = force_H + pressure * delta_Z;
					force_V = force_V + pressure * delta_X;
					point -= 1;
				}
			}

			//Add last segment
			Z2 = Z1;
			X2 = X1;
			bottom_pressure = top_pressure;
			Z1 = top_Z;
			X1 = top_X;
			top_pressure = 9.807 * (WS_Z - Z1);
			delta_X = X2 - X1;
			delta_Z = Z1 - Z2;
			pressure = 0.5 * (bottom_pressure + top_pressure);
			force_H = force_H + pressure * delta_Z;
			force_V = force_V + pressure * delta_X;

			//Compute confining force and angle
			confining_force[i] = sqrt(force_H * force_H + force_V * force_V);
			if (force_V == 0) {
				confining_angle[i] = 0.5 * M_PI;
			}
			else {
				confining_angle[i] = atan(force_H / force_V);
				if (confining_angle[i] < 0) { confining_angle[i] = confining_angle[i] + M_PI; }
			}

			//If top of soil layer is above the water surface elevation exit the for loop
			if (WS_Z <= bank_layer_Z[i]) { break; }
		}
	}
}

double compute_FoS(vector<vector<double>>& bank_geom, double& base_Z, double& base_X,
	double& failure_angle, vector<double>& unit_weight, vector<double>& cohesion, vector<double>& Phi,
	vector<double>& bank_layer_Z, vector<double>& bank_layer_X) {

	double FoS = 0;
	double Phib = 15 * M_PI / 180.0;

	//Find number of layers in failure plane block
	int num_layers = calc_num_layers(bank_geom, base_Z, bank_layer_Z);

	//Calculate intersect of failure plane and floodplain
	int FP_point = 21;
	vector<double> FP_intersect = calc_FP_shear_intersect(bank_geom, base_Z, base_X, failure_angle,
		FP_point);

	//Prevent small failures of top bank portion
	if (pow(pow(FP_intersect[0] - base_X, 2) + pow(FP_intersect[1] - base_Z, 2), 0.5) < 0.02) {
		FoS = 99999999;
		//GO TO END OF FUNCTION AND NEXT ITERATION
		return(FoS);
	}

	//Compute confining force and angle
	//Initialize confining force and angle vectors
	vector<double> confining_angle = { 0, 0 };
	vector<double> confining_force = { 0, 0 };
	calc_water_force(bank_geom, base_X, base_Z, confining_force, confining_angle, bank_layer_Z,
		bank_layer_X);
	//cout << "Driving Confining Force: " << confining_force[0] *
	//	sin(confining_angle[0] - failure_angle) << " " << confining_force[1] *
	//	sin(confining_angle[1] - failure_angle) << "\n";
	//cout << "Resisting Confining Force: " << confining_force[0] *
	//	cos(confining_angle[0] - failure_angle) * tan(Phi[0]) << " " <<
	//	confining_force[1] * cos(confining_angle[1] - failure_angle) * tan(Phi[1]) << "\n";

	//Compute the pore-water force of the layers
	//Create and adjusted unit weight vector to send to pore_water_force() and weight().
	vector<double> adj_unit_weight(2);
	adj_unit_weight[0] = unit_weight[0];
	adj_unit_weight[1] = unit_weight[1];
	vector<double> pore_water_force = calc_pore_water_force(base_Z, failure_angle, Phi,
		Phib, bank_geom, FP_intersect[1], bank_layer_Z, adj_unit_weight);
	//cout << "Pore force: " << pore_water_force[0] << " " << pore_water_force[1] << "\n";


	//Calculate weight of the layers
	vector<double> weight = calc_weight(bank_geom, FP_intersect, base_X, base_Z, adj_unit_weight,
		bank_layer_Z, bank_layer_X, failure_angle, FP_point);
	//cout << "weight: " << weight[0] << " " << weight[1] << "\n";

	//Intialize forces
	double sum_resisting_forces = 0, sum_driving_forces = 0;

	for (int i = 0; i < num_layers; i++) {
		double resisting_forces = (weight[i] * cos(failure_angle) + confining_force[i] *
			cos(confining_angle[i] - failure_angle)) * tan(Phi[i]) - pore_water_force[i];
		if (0 > resisting_forces) { resisting_forces = 0; }

		sum_resisting_forces += resisting_forces;
		sum_driving_forces += weight[i] * sin(failure_angle) - confining_force[i] *
			sin(confining_angle[i] - failure_angle);

		//Now add remaining resisting components (c'L + Pcos(a-b) tan f')

		//Adjust layer elevation if necessary
		double layer_Z = 0;
		if (i == 0) {
			layer_Z = max(FP_intersect[1], bank_layer_Z[i]);
		}
		else {
			layer_Z = bank_layer_Z[i];
		}

		//Compute length of failure plane within the layer
		double length = 0;
		if (FP_intersect[1] >= layer_Z) {
			length = (layer_Z - max(base_Z, bank_layer_Z[i + 1])) / sin(failure_angle);
		}
		else {
			length = (FP_intersect[1] - max(base_Z, bank_layer_Z[i + 1])) / sin(failure_angle);
		}
		if (0 > length) { length = 0; }

		sum_resisting_forces += cohesion[i] * length;
	}

	//If driving forces is zero, set FoS to large number
	if (0 >= sum_driving_forces) {
		FoS = 99999999;
	}
	else {
		FoS = sum_resisting_forces / sum_driving_forces;
	}

	if (FoS <= 0) { FoS = 99999999; }

	//cout << "Resisting: " << sum_resisting_forces << "\n";
	//cout << "Driving: " << sum_driving_forces << "\n";

	return(FoS);
}


int calc_num_layers(vector<vector<double>>& bank_geom, double& base_Z, vector<double>& bank_layer_Z) {

	//First find point on bank just above failure plane base
	//int point = 0;
	//while (bank_geom[1].size() > point) {
	//	if (bank_geom[1][point] < base_Z) { break; }
	//	point += 1;
	//}
	//point -= 1;

	//int num_layers = 0; //Number of layers in failure block
	//if (point > 16) {
	//	num_layers = 2; //toe and main bank
	//}
	//else {
	//	num_layers = 1; //just main bank
	//}
	int num_layers = 0;
	for (int i = 0; i < 2; ++i) {
		if (bank_layer_Z[i] > base_Z) { num_layers += 1; }
	}

	return(num_layers);
}

void calc_min_FoS(vector<vector<double>>& bank_geom, vector<double>& bank_layer_Z,
	vector<double>& bank_layer_X, vector<double>& cohesion, vector<double>& unit_weight,
	vector<double>& Phi, double& best_FoS, double& best_failure_angle,
	double& best_base_Z) {

	//Loop through 50 elevations for the failure plane base
	const int num_base_Z = 1; // changed to two to speed up calculations
	double new_FoS = 9999999, failure_angle = 0;

	////Find max and min elevations
	//double max_Z = -99999999;
	//double min_Z = 99999999;
	//for (int i = 0; i < (bank_geom[1].size()); ++i) {
	//	if (bank_geom[1][i] > max_Z) { max_Z = bank_geom[1][i]; }
	//	if (bank_geom[1][i] < min_Z) { min_Z = bank_geom[1][i]; }
	//}
	//best_base_Z = 0.5 * (max_Z + min_Z);
	////BSTEM trims max and min elevations by 1% of bank height but says this could be problematic...
	//double tmp1 = best_base_Z - 0.495 * (max_Z - min_Z);
	//double tmp2 = best_base_Z + 0.495 * (max_Z - min_Z);
	//min_Z = tmp1;
	//max_Z = tmp2;

	vector<double> new_Z = { bank_geom[1][16], bank_geom[1][21] };
	//Loop through the 50 iterations of possible base elevations and calculate FoS for each, storing
	//values if a smaller FoS is found
	for (int i = 0; i <= num_base_Z; ++i) {
		//double new_Z = min_Z + i * (max_Z - min_Z) / num_base_Z;
		bracket_brent(bank_geom, bank_layer_X, bank_layer_Z, cohesion, unit_weight, Phi,
			failure_angle, new_Z[i], new_FoS);
		if (best_FoS > new_FoS) {
			best_base_Z = new_Z[i];
			best_failure_angle = failure_angle;
			best_FoS = new_FoS;
		}
	}

	//Recalculate FoS??
	double base_X = set_bank_intersect_X(bank_geom, best_base_Z);
	best_FoS = compute_FoS(bank_geom, best_base_Z, base_X, best_failure_angle, unit_weight, cohesion,
		Phi, bank_layer_Z, bank_layer_X);
}

void bracket_brent(vector<vector<double>>& bank_geom, vector<double>& bank_layer_X,
	vector<double>& bank_layer_Z, vector<double>& cohesion, vector<double>& unit_weight,
	vector<double>& Phi, double& best_failure_angle, double& base_Z,
	double& best_FoS) {

	//Bracket angle and FoS values
	//Initialize angle and FoS arrays
	vector<double> angle(3);
	vector<double> FoS(3);

	//Find min angle
	angle[0] = calc_min_angle(base_Z, Phi, bank_geom, bank_layer_Z, bank_layer_X);

	//Calculate intersect of failure surface and bank face
	double base_X = set_bank_intersect_X(bank_geom, base_Z);

	//Calculate max angle
	angle[2] = calc_max_angle(base_Z, base_X, bank_geom);

	//If min angle > max angle, set FoS to big number and exit
	if (angle[0] > angle[2]) {
		best_failure_angle = angle[2];
		best_FoS = 99999999;
		return;
	}

	//Initialize angle 2 to mean of 1 and 3
	angle[1] = 0.5 * (angle[0] + angle[2]);

	//Comput initial FoS values
	for (int i = 0; i < 3; ++i) {
		FoS[i] = compute_FoS(bank_geom, base_Z, base_X, angle[i], unit_weight, cohesion, Phi,
			bank_layer_Z, bank_layer_X);
	}

	//Expand the interval as necessary to find new values with FoS[1] <= min(FoS[0], FoS[2])
	if (FoS[1] > FoS[0]) {
		if (FoS[2] >= FoS[0]) {
			//Expand interval to the left
			while (FoS[1] > FoS[0]) {
				if (0.5 * (angle[1] - angle[0]) < 8.0 * M_PI / 180.0) {
					//Cannot expand left: exit with the best values found so far
					best_failure_angle = angle[0];
					best_FoS = FoS[0];
					break;
				}
				//Expand left
				angle[2] = angle[1];
				FoS[2] = FoS[1];
				angle[1] = angle[0] + 0.5 * (angle[1] - angle[0]);
				//Comput FoS
				FoS[1] = compute_FoS(bank_geom, base_Z, base_X, angle[1], unit_weight, cohesion, Phi,
					bank_layer_Z, bank_layer_X);
			}
		}
		else {
			//We can't expand right: exit with the best values found so far
			best_failure_angle = angle[2];
			best_FoS = FoS[2];
		}
	}
	else if (FoS[1] > FoS[2]) {
		//Expand interval to the left
		while (FoS[1] > FoS[2]) {
			if (0.5 * (angle[2] - angle[1]) < 8.0 * M_PI / 180.0) {
				//Can't expand right: exit with best values found so far
				best_failure_angle = angle[2];
				best_FoS = FoS[2];
				break;
			}
			//Expand right
			angle[0] = angle[1];
			FoS[0] = FoS[1];
			angle[1] = angle[1] + 0.5 * (angle[2] - angle[1]);
			//Compute FoS
			FoS[1] = compute_FoS(bank_geom, base_Z, base_X, angle[1], unit_weight, cohesion, Phi,
				bank_layer_Z, bank_layer_X);
		}
		//We can't expand right: exit with the best values found so far
		best_failure_angle = angle[2];
		best_FoS = FoS[2];
	}
	else {
		best_failure_angle = angle[1];
		best_FoS = FoS[1];
	}

	/////////////////////////////////////////////////////////////////////
	//Brent algorithm
	double old_step = 0, step = 0, a = 0, b = 0, Q = 0, R = 0, P = 0, temp = 0, new_angle = 0,
		mean_angle = 0, new_FoS = 0;
	const int max_iterations = 1000;
	double tolerance = 0.25 * M_PI / 180.0; //tolerance of 0.25 degrees

	if (angle[0] < angle[2]) {
		a = angle[0];
		b = angle[2];
	}
	else {
		a = angle[2];
		b = angle[0];
	}

	//Re-intialize
	angle[1] = best_failure_angle;
	angle[0] = angle[2];
	FoS[1] = best_FoS;
	FoS[0] = FoS[1];

	//Main loop
	for (int i = 1; i <= max_iterations; ++i) {
		//Convergence check
		mean_angle = 0.5 * (a + b);
		if (abs(best_failure_angle - mean_angle) <= (2.0 * tolerance - 0.5 * (b - a))) { break; }

		//Construct trial parabolic fit
		if (abs(old_step) > tolerance) {
			Q = (best_failure_angle - angle[0]) * (best_FoS - FoS[1]);
			R = (best_failure_angle - angle[1]) * (best_FoS - FoS[0]);
			P = (best_failure_angle - angle[0]) * Q - (best_failure_angle - angle[1]) * R;
			Q = 2.0 * (Q - R);
			if (Q > 0) { P = -P; }
			Q = abs(Q);
			temp = old_step;
			old_step = step;

			//Test for acceptability of the parabolic fit. If it isn't successful, take the golden
			//section step into the large of the two segments
			if ((abs(P) >= abs(0.5 * Q * temp)) | (P <= Q * (a - best_failure_angle)) |
				(P >= Q * (b - best_failure_angle))) {
				if (best_failure_angle >= mean_angle) {
					old_step = a - best_failure_angle;
				}
				else {
					old_step = b - best_failure_angle;
				}
				step = 0.381966 * old_step;
			}
			else {
				//Take the parabolic step
				step = P / Q;
				new_angle = best_failure_angle + step;
				if ((new_angle - a < 2.0 * tolerance) | (b - new_angle < 2.0 * tolerance)) {
					double angle_diff = mean_angle - best_failure_angle;
					step = sign_fn(tolerance, angle_diff);
				}
			}
		}
		else {
			if (best_failure_angle >= mean_angle) {
				old_step = a - best_failure_angle;
			}
			else {
				old_step = b - best_failure_angle;
			}
			step = 0.381966 * old_step;
		}
		if (abs(step) >= tolerance) {
			new_angle = best_failure_angle + step;
		}
		else {
			new_angle = best_failure_angle + sign_fn(tolerance, step);
		}

		//Now calculate new FoS
		new_FoS = compute_FoS(bank_geom, base_Z, base_X, new_angle, unit_weight, cohesion, Phi,
			bank_layer_Z, bank_layer_X);

		//Assess quality of calculate FoS
		if (best_FoS >= new_FoS) {
			//Improvement
			if (new_angle >= best_failure_angle) {
				a = best_failure_angle;
			}
			else {
				b = best_failure_angle;
			}
			//Housekeeping
			shift_fn(angle[0], angle[1], best_failure_angle, new_angle);
			shift_fn(FoS[0], FoS[1], best_FoS, new_FoS);
		}
		else {
			if (new_angle < best_failure_angle) {
				a = new_angle;
			}
			else {
				b = new_angle;
			}
			if ((FoS[1] >= new_FoS) | (angle[1] == best_failure_angle)) {
				angle[0] = angle[1];
				angle[1] = new_angle;
				FoS[0] = FoS[1];
				FoS[1] = new_FoS;
			}
			else if ((FoS[0] >= new_FoS) | (angle[0] == best_failure_angle | angle[0] == angle[1])) {
				angle[0] = new_angle;
				FoS[0] = new_FoS;
			}
		}
	} //Next iterations
}

double sign_fn(double& x, double& y) {
	double Z = 0;
	if (y >= 0) {
		Z = abs(x);
	}
	else {
		Z = -abs(x);
	}
	return(Z);
}

void shift_fn(double& A, double& B, double& C, double& D) {
	A = B;
	B = C;
	C = D;
}