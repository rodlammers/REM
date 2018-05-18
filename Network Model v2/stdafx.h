// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"
#include "boost/multi_array.hpp"
#include "boost/math/tools/minima.hpp"
//#include "boost/filesystem.hpp"
#include "Incision Model.h"
#include "Transport Functions.h"
#include "Bank Erosion.h"
#include "Add_Functions.h"
#include "BSTEM Functions.h"
#include "Get_Inputs.h"
#include "Brent.h"

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <stdlib.h>
#include <numeric>
#include <cstring>
#include <omp.h>

extern vector<vector<string>> transport_err_names;
extern vector<vector<double>> transport_err_vals;

// TODO: reference additional headers your program requires here
