#pragma once

#include <cstdlib>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

struct bounds
{
	double left, right;
};

namespace function {
	typedef double (*function)(vector<double>);
}