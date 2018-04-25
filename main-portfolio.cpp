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

#define nDays 100
#define nSamples 100
#define alpha 0.4 // utility from pollution
#define beta 0.08 // eutrophic cost
#define reliability_threshold 0.85
#define inertia_threshold (-0.02)

double b, q, mu, sigma, delta, pcrit;

int nvars   = nDays;
int nobjs   = 4;
int nconsts = 1;

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
using namespace std;

ublas::vector<double> average_daily_P(nDays);
ublas::vector<double> discounted_benefit(nSamples);
ublas::vector<double> days_inertia_met(nSamples);
ublas::vector<double> days_pcrit_met(nSamples);

double modelmat[88][3] = {
		{6,18,362.3},
		{38,54,356},
		{73,75,209},
		{100,100,0},
		{14,5,781.5},
		{50,38,732},
		{83,74,118},
		{100,100,0},
		{24,0,76.8},
		{66,66,55},
		{99,80,3},
		{100,100,0},
		{0,9,2975},
		{19,13,2397},
		{28,24,2120},
		{30,30,0},
		{15,1,18438.7},
		{21,5,18033},
		{36,10,3608},
		{45,10,0},
		{6,6,445.5},
		{17,11,401},
		{40,25,297},
		{45,30,0},
		{20,7,8954.6},
		{47,34,8788},
		{93,80,7472},
		{100,100,0},
		{6,5,17.7},
		{26,50,14},
		{48,75,0},
		{50,80,0},
		{9,26,15.595},
		{36,38,14},
		{86,67,4},
		{100,80,0},
		{15,25,9.838},
		{21,28,7},
		{44,79,3},
		{60,80,0},
		{7,2,85.33},
		{15,6,58},
		{28,9,24},
		{30,10,0},
		{9,14,72.284},
		{25,27,71},
		{35,42,59},
		{50,50,0},
		{5,23,306.173},
		{16,59,239},
		{20,87,191},
		{30,100,0},
		{9,20,45.3},
		{44,41,39},
		{83,84,7},
		{100,100,0},
		{14,22,192.5},
		{27,33,188},
		{50,73,66},
		{50,80,0},
		{16,32,2415.8},
		{42,48,1810},
		{59,97,77},
		{80,100,0},
		{0,10,132.42},
		{5,41,112},
		{9,79,91},
		{10,100,0},
		{6,3,559.8},
		{13,13,393},
		{16,15,360},
		{20,20,0},
		{6,2,210.657},
		{11,45,146},
		{19,65,144},
		{20,75,0},
		{5,4,840.735},
		{24,23,746},
		{50,48,612},
		{50,50,0},
		{14,3,191.2},
		{50,28,177},
		{86,53,47},
		{100,60,0},
		{0,21,107.1},
		{0,55,90},
		{0,77,44},
		{0,100,0}
};

void lake_problem(double* vars, double* objs, double* consts) 
{
  zero(average_daily_P);
  zero(discounted_benefit);
  zero(days_inertia_met);
  zero(days_pcrit_met);

  for (int s = 0; s < nSamples; s++)
  {   
    // randomly generated natural phosphorous inflows
    ublas::vector<double> P_inflow = generateSOW(nDays, mu, sigma);

    double X = 0.0; // lake state 

    //implement the lake model from Carpenter et al. 1999
    for (int i = 0; i < nDays; i++)
    {
      // new state: previous state - decay + recycling + pollution
      X = X*(1-b) + pow(X,q)/(1+pow(X,q)) + vars[i] + P_inflow(i);

      average_daily_P(i) += X/nSamples;

      discounted_benefit(s) += alpha*vars[i]*pow(delta,i);

      if(i > 0 && vars[i]-vars[i-1] > inertia_threshold)
        days_inertia_met(s) += 1;

      if(X < pcrit)
        days_pcrit_met(s) += 1;

    }
  }

  // Calculate minimization objectives (defined in comments at beginning of file)
  objs[0] = vmax(average_daily_P);
  objs[1] = -1*vsum(discounted_benefit)/nSamples;
  objs[2] = -1*vsum(days_inertia_met)/((nDays-1)*nSamples);
  objs[3] = -1*vsum(days_pcrit_met)/(nDays*nSamples);

  consts[0] = max(0.0, reliability_threshold - (-1*objs[3]));

  average_daily_P.clear();
  discounted_benefit.clear();
  days_inertia_met.clear();
  days_pcrit_met.clear();
}

double root_function(double x) {
  return pow(x,q)/(1+pow(x,q)) - b*x;
}

bool root_termination(double min, double max) {
  return abs(max - min) <= 0.000001;
}

int main(int argc, char* argv[]) 
{  
  double vars[nvars];
  double objs[nobjs];
  double consts[nconsts];

  // initialize defaults
  b = 0.42;
  q = 2;
  mu = 0.02;
  sigma = 0.0017783;
  delta = 0.98;
  pcrit = -1;

  //read the command line arguments
  int opt;

  cout << modelmat[0][1]<<" "<< modelmat[0][2]<<" "<< modelmat[0][3]<<" "<<endl;

  while ((opt = getopt(argc, argv, "b:q:m:s:d:p:")) != -1) {
    switch (opt) {
      case 'b':
        b = atof(optarg);
        break;
      case 'q':
        q = atof(optarg);
        break;
      case 'm':
        mu = atof(optarg);
        break;
      case 's':
        sigma = atof(optarg);
        break;
      case 'd':
        delta = atof(optarg);
        break;
      case 'p':
        pcrit = atof(optarg);
        break;
      case '?':
      default:
        fprintf(stderr, "Unrecognized option\n");
        exit(EXIT_FAILURE);
    }
  }

  // compute the critical threshold if not supplied on command line
  if (pcrit < 0) {
    std::pair<double, double> result = tools::bisect(root_function, 0.01, 1.0, root_termination);
    pcrit = (result.first + result.second) / 2;
  }

  MOEA_Init(nobjs, nconsts);

  while (MOEA_Next_solution() == MOEA_SUCCESS) {
    MOEA_Read_doubles(nvars, vars);
    lake_problem(vars, objs, consts);
    MOEA_Write(objs, consts);
  }

  MOEA_Terminate();
  return EXIT_SUCCESS;
}
