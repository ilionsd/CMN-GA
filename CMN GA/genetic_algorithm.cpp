#pragma once

#include "genetic_algorithm.h"

vector<double>* decode(genetic_code* gcode, vector<double>* left_bounds, vector<double>* step_size, unsigned int dimension);

double euclidean_distance(vector<double> *v1, vector<double> *v2);

double select_rank(vector<double>* v, int rank, int size);
double select_median(vector<double>* v, int size);

vector<int>* nearest_neighbours(vector<double>* v, int n_min, int size);

int binary_search(vector<double> *v, int value);

genetic_algorithm::results* genetic_algorithm::cmn_ga(function::function func, vector<bounds>* bounds, parameters* parameters) {
	
	//-- общие параметны --
	int t_work;
	int t_start, t_finish;

	int generation_number = 0;
	//-- 3 знака после запятой --
	int precision = parameters->n_precision;

	int max_population = parameters->n_max_population;

	unsigned int dimension = bounds->size();
	vector<double> left_bound (dimension);
	vector<double> range_size (dimension);
	vector<double> step_size (dimension);
	double max_range = 0;
	for (int k = 0; k < dimension; ++k) {
		left_bound.at(k) = bounds->at(k).left;
		range_size.at(k) = bounds->at(k).right - bounds->at(k).left;
		if (max_range < range_size.at(k))
			max_range = range_size.at(k);
	}

	double stopping_criterion = pow(10, precision);
	chromosome_length chromosome_length = ceil(log(max_range * stopping_criterion * stopping_criterion));
	point_number point_count = 1 << chromosome_length;

	for (int k = 0; k < dimension; ++k) {
		step_size.at(k) = range_size.at(k) / (point_count - 1);
	}

	//--  --
	int local_optimims_final_count;
	vector<int>* local_optimums_final;
	vector<vector<double>*>* local_optimums_final_decoded;

	vector<double>* fitness0 = new vector<double> (max_population);
	vector<double>* fitness1;
	vector<double>* fitness2;
	vector<double>* fitness3 = new vector<double> (max_population);
	double min_fitness0 = 1e+100;
	double max_fitness0 = 0;

	vector<vector<double>*>* decoded = new vector<vector<double>*> (max_population);
	vector<vector<double>*>* invert_distance = new vector<vector<double>*> (max_population);

	t_start = clock();

#pragma region //-- Step 0: (Initialization) --
	cout << "Initialization...";
	vector<genetic_code*>* gcodes = new vector<genetic_code*> (max_population);
	int current_population = parameters->n_population;
	for (int k = 0; k < current_population; ++k) {
		gcodes->at(k) = genetic_code::generate(chromosome_length, dimension);
		decoded->at(k) = decode(gcodes->at(k), &left_bound, &step_size, dimension);
		fitness0->at(k) = func(*decoded->at(k));
		if (min_fitness0 > fitness0->at(k))
				min_fitness0 = fitness0->at(k);
	}

	for (int k = 0; k < current_population; ++k) 
		(*invert_distance)[k] = new vector<double> (max_population);

	for (int from = 0; from < current_population; ++from) {
		for (int to = from + 1; to < current_population; ++to) {
			double euc_distance = euclidean_distance((*decoded)[from], (*decoded)[to]);
			(*(*invert_distance)[from])[to] = 1.0 / euc_distance;
			(*(*invert_distance)[to])[from] = (*(*invert_distance)[from])[to];
		}
	}

	cout << " OK" << endl;

#pragma endregion

	double dot9_in_gpower = 1;

	bool population_overflow = false;

	while (true) {

#pragma region //-- Step 1: (Fitness Scaling) --

		++generation_number;
		cout << endl << "Population has " << current_population << " individuals" << endl;

		cout << "Fitness Scaling: generation " << generation_number << " ..." << endl; 

		int local_optimum_count = 0;
		vector<int>* local_optimums = new vector<int> (max_population);

		cout << "Local optimums detection ...";
		for (int k = 0; k < current_population; ++k) {
			vector<int>* neighbours = nearest_neighbours((*invert_distance)[k], parameters->n_min, current_population);
			bool is_local_optimum = true;
			for (vector<int>::iterator i = neighbours->begin(); is_local_optimum && i < neighbours->end(); ++i) {
				if ((*fitness0)[k] < (*fitness0)[*i])
					is_local_optimum = false;
			}
			if (is_local_optimum)
				(*local_optimums)[local_optimum_count++] = k; 
			delete neighbours;
		}
		cout << "Done. " << "Population has " << local_optimum_count << " optimums." << endl;
		cout << "Fitness 3 Calculating ... ";
		double fitness3_summ = 0;
		for (int ind = 0; ind < current_population; ++ind) {
			
			bool is_local_optimum = false;
			for (int k = 0; k < local_optimum_count && !is_local_optimum; ++k) 
				is_local_optimum = (ind == (*local_optimums)[k]);

			if (is_local_optimum) {
				//-- Если локальный оптимум, то масштабирование не нужно. Просто 1 --
				(*fitness3)[ind] = 1.0;
			}
			else {
				//-- Если не локальный оптимум, то масштабируем --
				//cout << "Fitness 1 Calculation ... ";
				fitness1 = new vector<double> (local_optimum_count);
				int local_optimums_fitter_count = 0;
				for (int k = 0; k < local_optimum_count; ++k) {
					if ((*fitness0)[ind] < (*fitness0)[(*local_optimums)[k]]) {
						double diff_fitness0 = (*fitness0)[(*local_optimums)[k]] - min_fitness0;
						(*fitness1)[k] = ((*fitness0)[ind] - min_fitness0) / diff_fitness0;
					}
					else {
						(*fitness1)[k] == 1.0;
					}
				}
				//cout << "Done" << endl;
				//cout << "Median Selection ... ";
				vector<double>* fitness_temp = new vector<double> (*fitness1);
				double median = select_median(fitness_temp, local_optimum_count);
				delete fitness_temp;
				//cout << "Done" << endl;
				double power = log(0.5) / log(median);
				//cout << "The Power is " << power << endl;
				//cout << "Fitness 2 Calculation ... ";
				fitness2 = new vector<double> (local_optimum_count);
				for (int k = 0; k < local_optimum_count; ++k) {
					(*fitness2)[k] = pow((*fitness1)[k], power);
				}
				//cout << "Done" << endl;
				
				double invert_distance_summ = 0;
				double invert_distance_weighted_summ = 0;
				for (int k = 0; k < local_optimum_count; ++k) {
					invert_distance_summ += (*(*invert_distance)[ind])[(*local_optimums)[k]];
					invert_distance_weighted_summ += (*(*invert_distance)[ind])[(*local_optimums)[k]] * (*fitness2)[k];
				}
				(*fitness3)[ind] = invert_distance_weighted_summ / invert_distance_summ;
				delete fitness1;
				delete fitness2;
			}
			fitness3_summ += (*fitness3)[ind];
		}
		cout << "Done" << endl;
		cout << " OK" << endl;

#pragma endregion

#pragma region //-- Step 4: (New generation) --

		cout << "Stopping Criterion: generation " << generation_number << " ..."; 

		//-- Stopping criterion --
		
		bool is_overall_converged = true;
		bool is_local_converged = false;
		for (int optimum_number = 0; optimum_number < local_optimum_count && is_overall_converged; ++optimum_number) {
			for (int k = 0; k < current_population && !is_local_converged; ++k) {
				if ((*(*invert_distance)[(*local_optimums)[optimum_number]])[k] > stopping_criterion)
					is_local_converged = true;
			}
			is_overall_converged = is_local_converged && is_overall_converged;
		}

		//-- Алгоритм сошёлся --
		if (is_overall_converged) {
			local_optimums_final = local_optimums;
			local_optimims_final_count = local_optimum_count;
			break;
		}

		//-- Алгоритм достиг максимальной популяции --
		if (population_overflow) {
			local_optimums_final = local_optimums;
			local_optimims_final_count = local_optimum_count;
			break;
		}

		delete local_optimums;

		cout << " OK" << endl;

#pragma endregion

		dot9_in_gpower *= 0.9;

#pragma region //-- Step 2: (Crossover) --

		cout << "Crossover: generation " << generation_number << " ... " << endl; 

		vector<int>* P1 = new vector<int> (parameters->n_crossover);
		vector<int>* P2 = new vector<int> (parameters->n_crossover);
		//-- P1: FPS --
		cout << "P1 Selecting: Probability Calculating ... ";
		vector<double>* probability = new vector<double> (current_population);
		(*probability)[0] = (*fitness3)[0] / fitness3_summ * RAND_MAX;
		for (int k = 1; k < current_population - 1; ++k) {
			(*probability)[k] = (*fitness3)[k] / fitness3_summ * RAND_MAX + (*probability)[k - 1];
		}
		(*probability)[current_population - 1] = RAND_MAX;
		cout << "P1 Selecting: Select ... ";
		for (int k = 0; k < parameters->n_crossover; ++k) {
			int rand_value = rand();
			int indIndex = binary_search(probability, rand_value);
			(*P1)[k] = indIndex;
		}
		delete probability;

		//-- P2: PPS --
		cout << "P2 Selecting ... ";
		for (int k = 0; k < parameters->n_crossover; ++k) {
			int max_proximitest = (double)rand() / (RAND_MAX + 1) * current_population;
			for (int i = 1; i < parameters->n_crowd; ++i) {
				int new_crowd = (double)rand() / (RAND_MAX + 1) * current_population; 
				if ((*(*invert_distance)[(*P1)[k]])[max_proximitest] < (*(*invert_distance)[(*P1)[k]])[new_crowd]) 
					max_proximitest = new_crowd;
			}
			(*P2)[k] = max_proximitest;
		}
		cout << "Selecting Done" << endl;
		//-- New Individuals --
		
		//int try_number = 0;
		cout << "Crossover ... " << endl;
		vector<bool>* p_used = new vector<bool>(parameters->n_crossover, false);
		int p_used_count = 0;
		cout << "Try ";
		for (int try_number = 0; try_number < parameters->n_try && !population_overflow && p_used_count < parameters->n_crossover; ++try_number) {
			cout << try_number << " ";
			for (int p_index = 0; p_index < parameters->n_crossover && try_number < parameters->n_try; ++p_index) {

				if ((*p_used)[p_index]) 
					continue;

				genetic_code *crossover_code = genetic_code::crossover((*gcodes)[(*P1)[p_index]], (*gcodes)[(*P2)[p_index]]);

				vector<double> *crossover_code_decoded = decode(crossover_code, &left_bound, &step_size, dimension);
				vector<double> *crossover_code_distance = new vector<double> (max_population);
				int min_distance_index = 0;
				for (int to = 0; to < current_population; ++to) {
					(*crossover_code_distance)[to] = euclidean_distance(crossover_code_decoded, (*decoded)[to]);
					if ((*crossover_code_distance)[min_distance_index] > (*crossover_code_distance)[to])
						min_distance_index = to;
				}
				double r_min = 0.08 * (1.001 - (*fitness3)[min_distance_index] * (1.0 - 0.5 * dot9_in_gpower));

				if (r_min < (*crossover_code_distance)[min_distance_index]) {
					//-- Добавляем --
					(*gcodes)[current_population] = crossover_code;
					(*decoded)[current_population] = crossover_code_decoded;
					(*fitness0)[current_population] = func(*crossover_code_decoded);
					if (min_fitness0 > (*fitness0)[current_population])
						min_fitness0 = (*fitness0)[current_population];
					for (int k = 0; k < crossover_code_distance->size(); ++k) {
						(*crossover_code_distance)[k] = 1 / (*crossover_code_distance)[k];
					}
					(*invert_distance)[current_population] = crossover_code_distance;
					for (int k = 0; k < current_population; ++k) 
						(*(*invert_distance)[k])[current_population] = (*crossover_code_distance)[k];
					++current_population;
					(*p_used)[p_index] = true;

					++p_used_count;

					if (current_population == max_population) {
						population_overflow = true;
						break;
					}
				}
				else {
					//-- Слишком близко --
					delete crossover_code_distance;
					delete crossover_code_decoded;
					delete crossover_code;
				}
			}
		}
		cout << endl;
		delete p_used;

		cout << " OK" << endl;

#pragma endregion


#pragma region //-- Step 3: (Mutation) --
		cout << "Mutation: generation " << generation_number << " ..." << endl;
		cout << "Try ";
		for (int try_number = 0; try_number < parameters->n_mutation && !population_overflow; ++try_number) {
			cout << try_number << " ";
			int mutant_index = (double)rand() / (RAND_MAX + 1) * current_population;
			genetic_code* mutant_code = genetic_code::mutate_normal((*gcodes)[mutant_index]);

			vector<double> *mutant_code_decoded = decode(mutant_code, &left_bound, &step_size, dimension);

			vector<double> *mutant_code_distance = new vector<double> (max_population);
			int min_distance_index = 0;
			for (int to = 0; to < current_population; ++to) {
				(*mutant_code_distance)[to] = euclidean_distance(mutant_code_decoded, (*decoded)[to]);
				if ((*mutant_code_distance)[min_distance_index] > (*mutant_code_distance)[to])
					min_distance_index = to;
			}
			double r_min = 0.08 * (1.001 - (*fitness3)[min_distance_index] * (1.0 - 0.5 * dot9_in_gpower));

			if (r_min < (*mutant_code_distance)[min_distance_index]) {
				//-- Добавляем --
				(*gcodes)[current_population] = mutant_code;
				(*decoded)[current_population] = mutant_code_decoded;
				(*fitness0)[current_population] = func(*mutant_code_decoded);
				if (min_fitness0 > (*fitness0)[current_population])
					min_fitness0 = (*fitness0)[current_population];
				for (int k = 0; k < mutant_code_distance->size(); ++k) {
					(*mutant_code_distance)[k] = 1 / (*mutant_code_distance)[k];
				}
				(*invert_distance)[current_population] = mutant_code_distance;
				for (int k = 0; k < current_population; ++k) 
					(*(*invert_distance)[k])[current_population] = (*mutant_code_distance)[k];
				++current_population;

				if (current_population == max_population) {
					population_overflow = true;
					break;
				}
			}
			else {
				//-- Слишком близко --
				delete mutant_code_distance;
				delete mutant_code_decoded;
				delete mutant_code;
			}
		}
		cout << " OK" << endl;
#pragma endregion

	}

	cout << " OK" << endl << "Finalazing ..."; 

	local_optimums_final_decoded = new vector<vector<double>*> (local_optimims_final_count);
	for (int k = 0; k < local_optimims_final_count; ++k) {
		(*local_optimums_final_decoded)[k] = new vector<double> ( *(*decoded)[(*local_optimums_final)[k]] );
	}

	//--  --
	t_finish = clock();

	t_work = (double)(t_finish - t_start) / CLOCKS_PER_SEC;
	for (int k = 0; k < current_population; ++k) {
		delete (*gcodes)[k];
		delete (*decoded)[k];
		delete (*invert_distance)[k];
	}
	delete gcodes;
	delete decoded;
	delete fitness0;
	delete fitness3;
	delete invert_distance;
	delete local_optimums_final;

	genetic_algorithm::results* results = new genetic_algorithm::results();
	results->optimums = local_optimums_final_decoded;
	results->work_time = t_work;
	results->final_population_size = current_population;
	results->final_generation_number = generation_number;

	cout << " OK" << endl;

	return results;
}

vector<double>* decode(genetic_code* gcode, vector<double>* left_bounds, vector<double>* step_size, unsigned int dimension) {
	vector<point_number>* numbers = genetic_code::decode(gcode);
	vector<double>* coordinates = new vector<double> (dimension);
	for (unsigned int k = 0; k < dimension; ++k) {
		coordinates->at(k) = left_bounds->at(k) + numbers->at(k) * step_size->at(k);
	}
	delete numbers;
	return coordinates;
}

double euclidean_distance(vector<double> *v1, vector<double> *v2) {
	int dimension = v1->size();
	double summ = 0;
	for (int k = 0; k < dimension; ++k) {
		double summand = (*v1)[k] - (*v2)[k];
		summ += summand * summand;
	}
	return sqrt(summ);
}

vector<int>* nearest_neighbours(vector<double>* v, int n_min, int size) {
	int rank = size - n_min - 1;

	vector<double>* temp_v = new vector<double> (size);
	for (int k = 0; k < size; ++k) 
		(*temp_v)[k] = (*v)[k];
	double rank_value = select_rank(temp_v, rank, size);
	delete temp_v;
	vector<int>* nearest = new vector<int> (n_min);
	int index = 0;
	for (int k = 0; k < size; ++k) {
		if ((*v)[k] > rank_value) {
			(*nearest)[index++] = k;
		}
	}
	return nearest;
}

#pragma region Quick Select

void swap(vector<double>* v, int i, int j) {
	double temp = v->at(i);
	(*v)[i] = (*v)[j];
	(*v)[j] = temp;
}

int partition(vector<double>* v, int left, int right, int pivotIndex) {
	double pivotValue = (*v)[pivotIndex];
	swap(v, pivotIndex, right);
	int storeIndex = left;
	for (int k = left; k < right; ++k) {
		if ((*v)[k] < pivotValue) {
			swap(v, storeIndex, k);
			++storeIndex;
		}
	}
	swap(v, right, storeIndex);
	return storeIndex;
}

double select_rank(vector<double>* v, int rank, int size) {
	int left = 0;
	int right = size - 1;

	while ( true ) {
		int pivotIndex = left + floor((double)rand() / (RAND_MAX + 1) * (right - left + 1));
		pivotIndex = partition(v, left, right, pivotIndex);
		if (rank == pivotIndex)
			return (*v)[rank];
		else if (rank < pivotIndex)
			right = pivotIndex - 1;
		else
			left = pivotIndex + 1;
	}
}

double select_median(vector<double>* v, int size) {
	int median = ceil(size / 2.0) - 1;
	return select_rank(v, median, size);

}

#pragma endregion

int binary_search(vector<double> *v, int value) {
	int left = 0;
	int right = v->size() - 1;
	
	while (left < right) {
		int middle = left + (right - left) / 2;
		if ((*v)[middle] < value) 
			left = middle + 1;
		else
			right = middle - 1;
	}
	return left;
}



