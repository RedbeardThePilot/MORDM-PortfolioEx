/* LakeProblem_5Obj_1Const_Stoch.cpp

 Riddhi Singh, May, 2014
 The Pennsylvania State University
 rus197@psu.edu

 Adapted by Tori Ward, July 2014
 Cornell University
 vlw27@cornell.edu

 Adapted by Jonathan Herman and David Hadka, Sept-Dec 2014
 Cornell University and The Pennsylvania State University

 A multi-objective represention of the lake model from Carpenter et al., 1999
 This simulation is designed for optimization with either Borg or the MOEAFramework.

 Stochasticity is introduced by natural phosphorous inflows.
 These follow a lognormal distribution with specified mu and sigma.

 Decision Vector
 vars : anthropogenic pollution flow at previous time step - Size 100, Bounds (0.0,0.1)

 Objectives
 1. minimize the maximum Phosphorous averaged over all lognormal draws in a given time period
 2. maximize expected benefit from pollution
 3. maximize the probability of meeting an inertia constraint
 4. maximize Reliability

 */

#include <stdio.h>
#include <unistd.h>
#include <sstream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/tools/roots.hpp>
#include "moeaframework.h"
#include "boostutil.h"
#include "generate-scenarios.h"

#include "modeldfn.h"

#define nPrograms 22
#define cost_threshold 30000000000

double bauScale, ssScale, costScale, budgetScale;

int nvars = nPrograms;
int nobjs = 3;
int nconsts = 1;

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
using namespace std;

void portfolio_problem(double* vars, double* objs, double* consts);

void portfolio_problem(double* vars, double* objs, double* consts) {
	double bau = 0, ss = 0, cost = 0;
	int optIdx; // Option Index
	vector<int> progOpts;
	progOpts.assign(vars, vars + nvars);

	for (int progIdx = 0; progIdx < nPrograms; progIdx++) {
		// new state: previous state - decay + recycling + pollution
		optIdx = progOpts[progIdx];

		bau += modelmat[4 * progIdx + optIdx][0];
		ss += modelmat[4 * progIdx + optIdx][1];
		cost += modelmat[4 * progIdx + optIdx][2];

	}

	// Calculate minimization objectives (defined in comments at beginning of file)
	objs[0] = bau * bauScale;
	objs[1] = ss * ssScale;
	objs[2] = cost * costScale;

	consts[0] = max(0.0, cost_threshold * budgetScale - cost);
}

//double root_function(double x) {
//	return pow(x, q) / (1 + pow(x, q)) - b * x;
//}

bool root_termination(double min, double max) {
	return abs(max - min) <= 0.000001;
}

int main(int argc, char* argv[]) {
	double vars[nvars];
	double objs[nobjs];
	double consts[nconsts];

	/* Initialize defaults */
	//Assume noble intent, polling of DMs is an accurate:
	bauScale = 1.0;
	ssScale = 1.0;
	costScale = 1.0;
	budgetScale = 1.0;

	/* read the command line arguments, if present (for openMORDM) */
	int opt;

	while ((opt = getopt(argc, argv, "b:q:m:s:d:p:")) != -1) {
		switch (opt) {
		case 'b': //Business as Usual Scale
			bauScale = atof(optarg);
			break;
		case 's':
			ssScale = atof(optarg);
			break;
		case 'c': //Cost Scale
			costScale = atof(optarg);
			break;
		case 'f': //Funds/Budget Scale
			costScale = atof(optarg);
			break;
		case '?':
		default:
			fprintf(stderr, "Unrecognized option\n");
			exit(EXIT_FAILURE);
		}
	}

	MOEA_Init(nobjs, nconsts);

	while (MOEA_Next_solution() == MOEA_SUCCESS) {
		MOEA_Read_doubles(nvars, vars);
		portfolio_problem(vars, objs, consts);
		MOEA_Write(objs, consts);
	}

	MOEA_Terminate();
	return EXIT_SUCCESS;
}
