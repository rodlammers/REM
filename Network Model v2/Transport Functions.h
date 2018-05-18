#ifndef TRANSPORT_FUNCTIONS_H
#define TRANSPORT_FUNCTIONS_H

#include "stdafx.h"
#include <string>
//using namespace std;

double qs_calc(const double& omega, const double& omega_c, const double& Ds, const double& q, const double& ps, 
	const string& type, const double& transport_factor);

#endif

