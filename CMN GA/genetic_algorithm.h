#pragma once

#include "function.h"
#include "genetic_code.h"

#include <vector>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <new>

using namespace std;
using namespace genetic;

namespace genetic_algorithm {
	//-- параметры алгоритма --
	struct parameters	{
		int n_population;
		int n_crossover;
		int n_mutation;
		int n_min;
		int n_crowd;
		int n_try;
		int n_max_population;
		int n_precision;
	};

	struct results {
		vector<vector<double>*>* optimums;
		unsigned int work_time;
		int final_population_size;
		int final_generation_number;
	};

	struct cresults {
		vector<double>* optimums;
		unsigned int work_time;
		int final_generation_number;
	};

	results* cmn_ga(function::function func, vector<bounds>* bounds, parameters* parameters);

	cresults* conventional_ga(function::function func,vector<bounds>* bounds, parameters* parameters);
}