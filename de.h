#pragma once
#include "global.h"
#include "bsptree.h"

struct DE
{
	vector<vector<double>> pop;
	vector<double> fitness;

	vector<vector<double>> trial;
	vector<double> f_trial;

	vector<double> best_solution;
	double best_fitness;

	double lb;
	double ub;
};

double GetRandom();

void InitializeDE(vector<DE> &de);

double Evaluate(vector<double> data, vector<int> &cur_evals, int t, int test_num);

void SetBound(vector<DE> &de, int test_num);

void InitializePop(DE &cur_de);

void EvaluatePop(BSPTree &cur_tree, vector<DE> &de, vector<vector<double>> &pop, vector<double> &fitness, vector<int> &cur_evals, int t, int test_num);

void Evolution(DE &cur_de, double f, double cr, int t);

void Selection(DE &cur_de);

void Transfer(vector<DE> &de, int t);

void LocalSearch(BSPTree &cur_tree, int idx, vector<DE> &de, int i_rand, int t, vector<int> cur_evals, int test_num);