#pragma once

#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "function.h"
#include "genetic_code.h"
#include "genetic_algorithm.h"

using namespace std;

double F0(vector<double> x);
double F1(vector<double> x);
double F2(vector<double> x);
double F4(vector<double> x);
double F5(vector<double> x);

int main(int argc, char *argv[]) {

	srand((unsigned) time( NULL ));

	//freopen("input.txt", "r", stdin);
	//freopen("log.txt", "w", stdout);

	ifstream input_file;
	ofstream output_file;
	input_file.open("input.txt");
	output_file.open("output.txt");

	/*
	int dimension, chromosome_length, codes_count;

	input_file >> dimension >> chromosome_length >> codes_count;

	for (int k = 0; k < codes_count; ++k) {
		genetic::genetic_code* gcode = genetic::genetic_code::generate(chromosome_length, dimension);
		for (vector<gene>::const_iterator i = gcode->gray_code->begin(); i < gcode->gray_code->end(); ++i)
			output_file << (int)(*i) << " ";
		output_file << "-> ";
		vector<point_number>* decoded = genetic::genetic_code::decode(gcode);
		delete gcode;
		for (vector<genetic::point_number>::const_iterator i = decoded->begin(); i < decoded->end(); ++i)
			output_file << (int)(*i) << "	";
		output_file << "->	";
		gcode = genetic::genetic_code::encode(chromosome_length, decoded);
		for (vector<gene>::const_iterator i = gcode->gray_code->begin(); i < gcode->gray_code->end(); ++i)
			output_file << (int)(*i) << " ";
		output_file << endl;
		delete gcode;
		delete decoded;
	}
	*/
	genetic_algorithm::parameters* params = new genetic_algorithm::parameters();
	input_file >> params->n_population;
	input_file >> params->n_crossover;
	input_file >> params->n_mutation;
	input_file >> params->n_min;
	input_file >> params->n_crowd;
	input_file >> params->n_try;
	input_file >> params->n_max_population;
	input_file >> params->n_precision;

	vector<bounds>* func_bounds = new vector<bounds> (2);
	(*func_bounds)[0].left = -40;
	(*func_bounds)[0].right = 40;

	(*func_bounds)[1].left = -40;
	(*func_bounds)[1].right = 40;


	function::function objective_function = &F4;

	genetic_algorithm::results* results = genetic_algorithm::cmn_ga(objective_function, func_bounds, params);

	output_file << "Завершено за " << results->work_time << "c. с размером популяции " << results->final_population_size << " и " << results->final_generation_number << " поколения" << endl;
	for (int k = 0; k < results->optimums->size(); ++k) {
		double value = objective_function(*(*results->optimums)[k]);
		for (int h = 0; h < (*results->optimums)[k]->size(); ++h) {
			output_file << (*(*results->optimums)[k])[h] << " ";
		}
		output_file << "-> " << value << endl;
	}

	delete params;
	delete func_bounds;
	delete results;
	
	input_file.close();
	output_file.close();

	return 0;
}


double F1(vector<double> x) {
	double arg = 5.1 * M_PI * x[0] + 0.5;
	return pow(sin(arg), 6);
}

double F0(vector<double> x) {
	return sin(M_PI * x[0]) + 0.1 * sin(100 * M_PI * x[0]);
}

double F2(vector<double> x) {
	double sin_arg = 5.1 * M_PI * x[0] + 0.5;
	double exp_arg = -(4*log(2)*pow(x[0] - 0.0667, 2)) / 0.64;
	return exp(exp_arg) * pow(sin(sin_arg), 6);
}

double F4(vector<double> x) {
	double A[5] = {-20, -5, 0, 30, 30};
	double B[5] = {-20, -25, 30, 0, -30};
	double H[5] = {0.4, 0.2, 0.7, 1.0, 0.05};
	double W[5] = {0.02, 0.5, 0.01, 2.0, 0.1};
	double summ = 0;
	for (int k = 0; k < 5; ++k) {
		summ += H[k] / (1 + W[k] * (pow(x[0] - A[k], 2) + pow(x[1] - B[k], 2)));
	}
	return summ;
}

double F5(vector<double> x) {
	return sin(M_PI * x[0] / 100.0);
}