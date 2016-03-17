#include "genetic_code.h"

genetic::genetic_code::genetic_code() {
}

genetic::genetic_code::genetic_code(genetic::chromosome_length chromosome_length, vector<genetic::gene> *genetic_code) {
	this->chromosomeLength = chromosome_length;
	this->gray_code = new vector<gene> (*genetic_code);
}

genetic::genetic_code::~genetic_code() {
	delete (this->gray_code);
}

void genetic::genetic_code::init_rand() {
	srand((unsigned) time( NULL ));
}

genetic::genetic_code* genetic::genetic_code::generate(chromosome_length chromosome_length, int dimension) {
	vector<gene> *gcode = new vector<gene> (chromosome_length * dimension);
	for (vector<gene>::iterator k = gcode->begin(); k < gcode->end(); ++k) 
		(*k) = floor((double)rand() / RAND_MAX + 0.5);
	return new genetic_code(chromosome_length, gcode);
}

genetic::genetic_code* genetic::genetic_code::encode(chromosome_length chromosome_length, vector<point_number>* numbers) {
	unsigned int dimension = numbers->size();
	unsigned int code_size = chromosome_length * dimension; 

	vector<gene>* gray_codes = new vector<gene> (code_size);

	for (unsigned int k = 0; k < dimension; ++k) {
		point_number gray_code = (*numbers)[k] ^ ((*numbers)[k] >> 1);
		unsigned int chromosome_ptr = k * chromosome_length;
		for (unsigned int gene_ptr = chromosome_length - 1; gene_ptr > 0 ; --gene_ptr) {
			(*gray_codes)[chromosome_ptr + gene_ptr] = (gene)(gray_code & 0x1);
			gray_code = gray_code >> 1;
		}
		(*gray_codes)[chromosome_ptr] = (gene)(gray_code & 0x1);
	}
	return new genetic_code(chromosome_length, gray_codes);
}

vector<genetic::point_number>* genetic::genetic_code::decode(genetic_code *gcode) {
	unsigned int code_size = gcode->gray_code->size();
	chromosome_length chromosome_length = gcode->chromosomeLength;
	unsigned int dimension = code_size / chromosome_length;

	vector<point_number> *numbers = new vector<point_number> (dimension, 0);
	vector<gene> binary (chromosome_length, 0);

	for (unsigned int k = 0; k < dimension; ++k) {
		point_number number = 0;
		unsigned int chromosome_ptr = k * chromosome_length;

		for (unsigned int digit_ptr = 0; digit_ptr < chromosome_length; ++digit_ptr) {
			for (unsigned int shift = digit_ptr; shift < chromosome_length; ++shift) {
				binary[shift] ^= (gcode->gray_code)->at(chromosome_ptr + shift - digit_ptr);
			}
			number <<= 1;
			number += binary.at(digit_ptr);
		}
		numbers->at(k) = number;
		binary.assign(chromosome_length, 0);
	}
	return numbers;
}

vector<double>* genetic::genetic_code::decode(genetic_code *gcode, vector<bounds> *bounds) {
	vector<point_number>* numbers = genetic_code::decode(gcode);
	unsigned int dimension = numbers->size();
	int point_count = 1 << gcode->chromosomeLength;
	vector<double>* coordinates = new vector<double> (dimension);
	for (unsigned int k = 0; k < dimension; ++k) {
		coordinates->at(k) = bounds->at(k).left + numbers->at(k) * (bounds->at(k).right - bounds->at(k).left) / point_count;
	}
	delete [] numbers;
	return coordinates;
}

genetic::genetic_code* genetic::genetic_code::crossover(genetic_code* parent_1, genetic_code* parent_2) {
	//-- Размеры кодов должны совпадать!!! ---
	unsigned int code_size = parent_1->gray_code->size(); 
	vector<gene>* new_code = new vector<gene> (code_size);
	for (int k = 0; k < code_size; ++k) {
		if (floor((double)rand() / RAND_MAX + 0.5)) 
			(*new_code)[k] = (*parent_1->gray_code)[k];
		else 
			(*new_code)[k] = (*parent_2->gray_code)[k];
	}
	genetic_code* offspring = new genetic_code(parent_1->chromosomeLength, new_code);
	delete new_code;
	return offspring;
}

double genetic::normal_rand(double expectation_value, double standart_deviation) {
	double u, v, s;
	do {
		u = (double)rand() / (RAND_MAX + 1) * 2 - 1;
		v = (double)rand() / (RAND_MAX + 1) * 2 - 1;
		s = u * u + v * v;
	} while (s > 1 || fabs(s) < 1e-9); 
	double gaussian_destribution = u * sqrt((-2.0) * log(s) / s);
	return standart_deviation * gaussian_destribution + expectation_value;
}

genetic::genetic_code* genetic::genetic_code::mutate(genetic_code* mutant) {
	int code_size = mutant->gray_code->size();
	vector<gene> new_code (code_size);
	for (int k = 0; k < code_size; ++k) {
		if ((double)normal_rand() / RAND_MAX < 0.5)
			new_code[k] = (*mutant->gray_code)[k];
		else 
			new_code[k] = !(*mutant->gray_code)[k];
	}
	return new genetic_code(mutant->chromosomeLength, &new_code);
}

genetic::genetic_code* genetic::genetic_code::mutate_normal(genetic_code* mutant) {
	int dimension = mutant->gray_code->size() / mutant->chromosomeLength;
	int dimendions_size = 1 << mutant->chromosomeLength;
	//-- 40% design space dimention--
	double standart_deviation = 0.4 * dimendions_size;
	vector<point_number>* decoded = genetic_code::decode(mutant);
	vector<point_number>* mutate = new vector<point_number> (dimension);
	for (int k = 0; k < dimension; ++k) {
		double real_value = genetic::normal_rand((double)(*decoded)[k], standart_deviation);
		(*mutate)[k] = floor(real_value + 0.5);
	}
	genetic_code* encoded = genetic_code::encode(mutant->chromosomeLength, mutate);
	delete decoded;
	delete mutate;
	return encoded;
}