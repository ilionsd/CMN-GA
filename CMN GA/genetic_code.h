#pragma once

#include "function.h"

#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

namespace genetic {
	typedef unsigned __int32 chromosome_length;
	typedef __int32 point_number;
	typedef __int8 gene;

	double normal_rand(double expectation_value = 0, double standart_deviation = 1);

	class genetic_code {
	public:
		vector<gene> *gray_code;
		chromosome_length chromosomeLength;

	private: 
		genetic_code(chromosome_length chromosome_length, vector<gene> *genetic_code);

	public:
		genetic_code();
		genetic_code(int dimension, double *coordinates);
		~genetic_code();
		static void init_rand();
		static genetic_code* generate(chromosome_length chromosome_length, int dimension);

		static genetic_code* encode(chromosome_length chromosome_length, vector<point_number>* numbers);

		static vector<point_number>* decode(genetic_code *gcode);
		static vector<double>* decode(genetic_code *gcode, vector<bounds> *bounds);

		static double euclidean_distance(genetic_code *gcode1, genetic_code *gcode2);

		static genetic_code* crossover(genetic_code* parent_1, genetic_code* parent_2);
		static genetic_code* mutate(genetic_code* mutant);
		static genetic_code* mutate_normal(genetic_code* mutant);
	};
}