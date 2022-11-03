#include "de.h"
#include "benchmark.h"
#include "bsptree.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <random>

double GetRandom()
{
	std::random_device seed;
	std::mt19937_64 engine(seed());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	return dist(engine);
	// return (double)rand() / (RAND_MAX + 1.0);
}

void InitializeDE(vector<DE> &de)
{
	for (int t = 0; t < de.size(); t++)
	{
		de[t].pop = vector<vector<double>>(POP_SIZE, vector<double>(DIMENSION));
		de[t].fitness = vector<double>(POP_SIZE, DBL_MAX);

		de[t].trial = vector<vector<double>>(POP_SIZE, vector<double>(DIMENSION));
		de[t].f_trial = vector<double>(POP_SIZE, DBL_MAX);

		de[t].best_solution = vector<double>(DIMENSION);
		de[t].best_fitness = DBL_MAX;
	}
}

double Evaluate(vector<double> data, vector<int> &cur_evals, int t, int test_num)
{
	cur_evals[t]++;
	double result = 0.0;
	switch (test_num)
	{
	case 1:
		if (t == 0)
			result = Griewank(data);
		else
			result = Rastrigin(data);
		break;
	case 2:
		if (t == 0)
			result = Ackley(data);
		else
			result = Rastrigin(data);
		break;
	case 3:
		if (t == 0)
			result = Ackley(data);
		else
			result = Schwefel(data);
		break;
	case 4:
		if (t == 0)
			result = Rastrigin(data);
		else
			result = Sphere(data);
		break;
	case 5:
		if (t == 0)
			result = Ackley(data);
		else
			result = Rosenbrock(data);
		break;
	case 6:
		if (t == 0)
			result = Ackley(data);
		else
			result = Weierstrass(data);
		break;
	case 7:
		if (t == 0)
			result = Rosenbrock(data);
		else
			result = Rastrigin(data);
	case 8:
		if (t == 0)
			result = Griewank(data);
		else
			result = Weierstrass(data);
		break;
	case 9:
		if (t == 0)
			result = Rastrigin(data);
		else
			result = Schwefel(data);
		break;
	default:
		printf("error problem number!");
		break;
	}

	return result;
}

void SetBound(vector<DE> &de, int test_num)
{
	switch (test_num)
	{

	case 1: // complete intersection with high similarity, Griewankand Rastrigin
		de[0].lb = -100;
		de[0].ub = 100;

		de[1].lb = -50;
		de[1].ub = 50;
		break;

	case 2: // complete intersection with medium similarity, Ackleyand Rastrigin
		de[0].lb = -50;
		de[0].ub = 50;

		de[1].lb = -50;
		de[1].ub = 50;
		break;

	case 3: // complete intersection with low similarity, Ackleyand Schwefel
		de[0].lb = -50;
		de[0].ub = 50;

		de[1].lb = -50;
		de[1].ub = 50;
		break;

	case 4: // partially intersection with high similarity, Rastriginand Sphere
		de[0].lb = -50;
		de[0].ub = 50;

		de[1].lb = -100;
		de[1].ub = 100;
		break;

	case 5: // partially intersection with medium similarity, Ackleyand Rosenbrock
		de[0].lb = -50;
		de[0].ub = 50;

		de[1].lb = -50;
		de[1].ub = 50;
		break;

	case 6: // partially intersection with low similarity, Ackley and Weierstrass
		de[0].lb = -50;
		de[0].ub = 50;

		de[1].lb = -0.5;
		de[1].ub = 0.5;
		break;

	case 7: // no intersection with high similarity, Rosenbrockand Rastrigin
		de[0].lb = -50;
		de[0].ub = 50;

		de[1].lb = -50;
		de[1].ub = 50;
		break;

	case 8: //% no intersection with medium similarity, Griewankand Weierstrass
		de[0].lb = -100;
		de[0].ub = 100;

		de[1].lb = -0.5;
		de[1].ub = 0.5;
		break;

	case 9: // % no overlap with low similarity, Rastriginand Schwefel
		de[0].lb = -50;
		de[0].ub = 50;

		de[1].lb = -500;
		de[1].ub = 500;
		break;

	default:
		break;
	}
}

void InitializePop(DE &cur_de)
{
	for (int i = 0; i < POP_SIZE; ++i)
		for (int j = 0; j < DIMENSION; ++j)
			cur_de.pop[i][j] = cur_de.lb + GetRandom() * (cur_de.ub - cur_de.lb);
	cur_de.best_fitness = DBL_MAX;
}

void EvaluatePop(BSPTree &cur_tree, vector<DE> &de, vector<vector<double>> &pop, vector<double> &fitness, vector<int> &cur_evals, int t, int test_num)
{
	for (int i = 0; i < POP_SIZE; ++i)
	{
		int idx = SearchNode(cur_tree, pop[i], de[t].lb, de[t].ub);
		if (idx == -1)
		{
			de[t].fitness[i] = Evaluate(pop[i], cur_evals, t, test_num);
			InsertNode(cur_tree, idx, pop[i], fitness[i]);
		}
		else
		{
			double dis = CalculateDistance(pop[i], cur_tree.node[i]);
			if (dis == 0)
			{
				LocalSearch(cur_tree, idx, de, i, t, cur_evals, test_num);
			}
			else
			{
				fitness[i] = Evaluate(pop[i], cur_evals, t, test_num);
				InsertNode(cur_tree, idx, pop[i], fitness[i]);
			}
		}
	}
}

void Evolution(DE &cur_de, double f, double cr, int t)
{
	for (int i = 0; i < POP_SIZE; ++i)
	{
		int r1, r2, r3;
		do
		{
			r1 = rand() % POP_SIZE;
			r2 = rand() % POP_SIZE;
			r3 = rand() % POP_SIZE;
		} while (r1 == i || r2 == i || r3 == i || r1 == r2 || r2 == r3);

		int r = rand() % DIMENSION;

		for (int j = 0; j < DIMENSION; ++j)
		{
			if (j == r || GetRandom() < cr)
			{
				cur_de.trial[i][j] = cur_de.pop[r1][j] + f * (cur_de.pop[r2][j] - cur_de.pop[r3][j]);
				if (cur_de.trial[i][j] > cur_de.ub || cur_de.trial[i][j] < cur_de.lb)
					cur_de.trial[i][j] = cur_de.lb + GetRandom() * (cur_de.ub - cur_de.lb);
			}
			else
			{
				cur_de.trial[i][j] = cur_de.pop[i][j];
			}
		}
	}
}

void Selection(DE &cur_de)
{
	for (int i = 0; i < POP_SIZE; ++i)
	{
		if (cur_de.f_trial[i] < cur_de.fitness[i])
		{
			cur_de.pop[i] = cur_de.trial[i];
			cur_de.fitness[i] = cur_de.f_trial[i];
		}

		if (cur_de.fitness[i] < cur_de.best_fitness)
		{
			cur_de.best_solution = cur_de.pop[i];
			cur_de.best_fitness = cur_de.fitness[i];
		}
	}
}

void Transfer(vector<DE> &de, int t)
{
	int a_t = (t + 1) % T;
	int r = rand() % POP_SIZE;
	de[t].pop[r] = de[a_t].best_solution;
}

void LocalSearch(BSPTree &cur_tree, int idx, vector<DE> &de, int i_rand, int t, vector<int> cur_evals, int test_num)
{
	int dim = de[0].pop[0].size();
	vector<double> data(dim);
	int l_idx = cur_tree.l_idx[cur_tree.p_idx[idx]];
	int r_idx = cur_tree.r_idx[cur_tree.p_idx[idx]];
	for (int i = 0; i < dim; i++)
	{
		// data[i] = cur_tree.sub_h[0][i] + GetRandom() * (cur_tree.sub_h[1][i] - cur_tree.sub_h[0][i]);
		data[i] = cur_tree.node[l_idx][i] + GetRandom() * (cur_tree.node[r_idx][i] - cur_tree.node[l_idx][i]);
	}
	int s = SearchNode(cur_tree, data, de[t].lb, de[t].ub);
	if (CalculateDistance(cur_tree.node[s], data))
	{
		double fitness = Evaluate(data, cur_evals, t, test_num);
		InsertNode(cur_tree, s, data, fitness);
		de[t].pop[i_rand] = data;
		de[t].fitness[i_rand] = fitness;
	}
}