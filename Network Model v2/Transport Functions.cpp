
#include "stdafx.h"
//#include "Transport Functions.h"
//#include <string>
using namespace std;

//qs_calc calculates volumetric unit sediment transport rates (m^2/s) using a variety of transport
//equations. These include new empirical, stream power based bedload and total load equations 
//(see Lammers and Bledsoe, 2017, ESPL). Two other bedload equations are also included 
//(Eaton and Church 2011, ESPL and Parker et al. 2012, River Res. Applic.) The user 
//specifies which transport equation to use. This choice is stored in the 'type' variable. Options 
//are "bedload", "total", "EC", "Parker", and "both". "Both" uses both the bedload and total load equations
//depending on the grain size being modeled (total load for <= 4 mm, and bedload for > 4 mm).

//INPUTS:
//omega: specific stream power (W/m^2)
//omega_c: critical specific stream power (W/m^2)
//Ds: grain size (m)
//q: unit discharge (m^3/m/s)
//ps: fraction of the bed material of the given grain size (dimensionless; 0-1)
//type: sediment transport equation to use

//OUTPUT:
//qs: volumetric unit sediment transport rate (m^3/m/s)

double qs_calc(const double& omega, const double& omega_c, const double& Ds, const double& q,
	const double& ps, const string& type, const double& transport_factor) {
	double qs = 0; //calculated volumetric unit sediment transport rate (m^3/m/s)
	double Qppm, omega_star, E_star;
	double omegac_star = 0.1; //constant critical dimensionless specific stream power
	double sg = 2.65; //specific gravity of sediment

	//Adjust sg based on percentage dacite and pumice - Used for North Fork Toutle River modeling
	//sg = 2.36 * pow(Ds * 1000, -0.036) * 0.3 + 1.44 * pow(Ds * 1000, -0.09) * 0.3 + 2.65 * 0.4;

	if (type == "bedload" | type == "total" | type == "both") {
		//Using bedload and/or total load equations. Transport only predicted if specific stream power
		//exceeds critical stream power for this grain size.
		if (omega > omega_c) {
			if (type == "bedload") {
				//Bedload equation - converts rate from kg/m/s to m^3/m/s
				qs = (omega - omega_c) * sqrt(omega - omega_c) / sqrt(Ds) / sqrt(q) * 1.43e-4 / (sg * 1000) * 
					ps * transport_factor;

				//Checks values for bedload equation against ranges used to fit the equation
				double qs2 = qs * (sg * 1000) / ps;
				transport_eq_check(omega, q, Ds, qs2, type);
			}
			else if (type == "total") {
				//Total load equation - PPM
				Qppm = (omega - omega_c) * sqrt(omega - omega_c) / Ds * pow(q, -5.0 / 6.0) * 0.0214 * transport_factor;
				
				if (Qppm > 1600000.0) {
					Qppm = 1600000.0;
				}

				//Total load in PPM is converted to m^3/m/s
				qs = Qppm * sg * 1e6 / (sg * 1e6 - (sg - 1) * Qppm) * q / (1000 * 1000 * sg) * ps;
				
				//Checks values for total load equation against ranges used to fit the equation
				transport_eq_check(omega, q, Ds, Qppm, "total");
			}
			else if (type == "both") {
				//Uses both total load and bedload equations depending on sediment grain size
				if (Ds <= 0.004) {
					Qppm = (omega - omega_c) * sqrt(omega - omega_c) / Ds * pow(q, -5.0 / 6.0) * 0.0214 * transport_factor;

					if (Qppm > 1600000.0) {
						Qppm = 1600000.0;
					}

					//Total load in PPM is converted to m^3/m/s
					qs = Qppm * sg * 1e6 / (sg * 1e6 - (sg - 1) * Qppm) * q / (1000 * 1000 * sg) * ps;
					
					//Checks values for total load equation against ranges used to fit the equation
					transport_eq_check(omega, q, Ds, Qppm, "total");
				}
				else {
					//Bedload equation - converts rate from kg/m/s to m^3/m/s
					qs = (omega - omega_c) * sqrt(omega - omega_c) / sqrt(Ds) / sqrt(q) * 1.43e-4 / (sg * 1000) * 
						ps * transport_factor;

					//Checks values for bedload equation against ranges used to fit the equation
					double qs2 = qs * (sg * 1000) / ps;
					transport_eq_check(omega, q, Ds, qs2, "bedload");
				}
			}
		}
	}
	else if (type == "EC") {
		//Eaton and Church (2011) equation
		if (omega == 0) {
			qs = 0;
		}
		else {
			omega_star = omega / (1000 * (sg - 1) * 9.81 * Ds * sqrt((sg - 1) * 9.81 * Ds));
			E_star = pow(0.92 - 0.25 * sqrt(omegac_star / omega_star), 9.0);
			qs = E_star * (omega / 9810.0) / (sg - 1) * ps * transport_factor;
		}
	}
	else {
		//Parker et al. (2012) equation
		if (omega == 0) {
			qs = 0;
		}
		else{
			double qb_star = 0;
			omega_star = omega / (1000 * (sg - 1) * 9.81 * Ds * sqrt((sg - 1) * 9.81 * Ds));
			if (omega_star < 0.25) {
				qb_star = 100 * pow(omega_star, 6);
			}
			else {
				qb_star = 0.2 * omega_star * sqrt(omega_star);
			}
			qs = qb_star * (1000 * 1.65 * 9.81 * Ds * sqrt(1.65 * 9.81 * Ds)) * transport_factor; //qs in N/m/s submerged weight
			qs = qs / 9.81 * sg / (sg - 1) / (sg * 1000) * ps; //Convert to m^3/m/s dry mass
		}
	}

	//Catch case where qs is really really small - used to prevent divide by zero errors
	if (qs < 1e-50) {
		qs = 0;
	}

	return qs;
}
