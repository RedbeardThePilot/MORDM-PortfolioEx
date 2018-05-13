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
#define cost_threshold 35000

double bauScale, ssScale, costScale, budgetScale;

int nvars = nPrograms;
int nobjs = 3;
int nconsts = 1;

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
using namespace std;

void portfolio_problem(double* vars, double* objs, double* consts, double* uncertainty);

void portfolio_problem(double* vars, double* objs, double* consts, double* uncertainty) {
	double bau = 0, ss = 0, cost = 0;
	int optIdx; // Option Index
	vector<int> progOpts;
	progOpts.assign(vars, vars + nvars);

	for (int progIdx = 0; progIdx < nPrograms; progIdx++) {
		optIdx = progOpts[progIdx];

		bau += uncertainty[progIdx] * modelmat[4 * progIdx + optIdx][0];
		ss += uncertainty[progIdx] * modelmat[4 * progIdx + optIdx][1];
		cost += uncertainty[progIdx] * modelmat[4 * progIdx + optIdx][2];

	}

	// Calculate minimization objectives (defined in comments at beginning of file)
	objs[0] = bau * bauScale;
	objs[1] = ss * ssScale;
	objs[2] = cost * costScale;

	consts[0] = max(0.0, cost - cost_threshold * budgetScale);
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
	double uncertainty[nPrograms];

	/* Initialize defaults */
	//Assume noble intent, polling of DMs is an accurate:
	bauScale = 1.0;
	ssScale = 1.0;
	costScale = 1.0;
	budgetScale = 1.0;
	fill_n(uncertainty, 22, 1);

	/* read the command line arguments, if present (for openMORDM) */
	int opt;

	while ((opt = getopt(argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:")) != -1) {
		switch (opt) {
		case 'w': //Business as Usual Scale
			bauScale = atof(optarg);
			break;
		case 'x':
			ssScale = atof(optarg);
			break;
		case 'y': //Cost Scale
			costScale = atof(optarg);
			break;
		case 'z': //Funds/Budget Scale
			costScale = atof(optarg);
			break;
		case 'a': //Uncertainty Multiplier for Program 1
			uncertainty[0] = atof(optarg);
			break;
		case 'b': //Uncertainty Multiplier for Program 2
			uncertainty[1] = atof(optarg);
			break;
		case 'c': //Uncertainty Multiplier for Program 3
			uncertainty[2] = atof(optarg);
			break;
		case 'd': //Uncertainty Multiplier for Program 4
			uncertainty[3] = atof(optarg);
			break;
		case 'e': //Uncertainty Multiplier for Program 5
			uncertainty[4] = atof(optarg);
			break;
		case 'f': //Uncertainty Multiplier for Program 6
			uncertainty[5] = atof(optarg);
			break;
		case 'g': //Uncertainty Multiplier for Program 7
			uncertainty[6] = atof(optarg);
			break;
		case 'h': //Uncertainty Multiplier for Program 8
			uncertainty[7] = atof(optarg);
			break;
		case 'i': //Uncertainty Multiplier for Program 9
			uncertainty[8] = atof(optarg);
			break;
		case 'j': //Uncertainty Multiplier for Program 10
			uncertainty[9] = atof(optarg);
			break;
		case 'k': //Uncertainty Multiplier for Program 11
			uncertainty[10] = atof(optarg);
			break;
		case 'l': //Uncertainty Multiplier for Program 12
			uncertainty[11] = atof(optarg);
			break;
		case 'm': //Uncertainty Multiplier for Program 13
			uncertainty[12] = atof(optarg);
			break;
		case 'n': //Uncertainty Multiplier for Program 14
			uncertainty[13] = atof(optarg);
			break;
		case 'o': //Uncertainty Multiplier for Program 15
			uncertainty[14] = atof(optarg);
			break;
		case 'p': //Uncertainty Multiplier for Program 16
			uncertainty[15] = atof(optarg);
			break;
		case 'q': //Uncertainty Multiplier for Program 17
			uncertainty[16] = atof(optarg);
			break;
		case 'r': //Uncertainty Multiplier for Program 18
			uncertainty[17] = atof(optarg);
			break;
		case 's': //Uncertainty Multiplier for Program 19
			uncertainty[18] = atof(optarg);
			break;
		case 't': //Uncertainty Multiplier for Program 20
			uncertainty[19] = atof(optarg);
			break;
		case 'u': //Uncertainty Multiplier for Program 21
			uncertainty[20] = atof(optarg);
			break;
		case 'v': //Uncertainty Multiplier for Program 22
			uncertainty[21] = atof(optarg);
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
		portfolio_problem(vars, objs, consts, uncertainty);
		MOEA_Write(objs, consts);
	}

	MOEA_Terminate();
	return EXIT_SUCCESS;
}
