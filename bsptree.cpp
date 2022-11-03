#include "bsptree.h"
#include "de.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <algorithm>

void InitializeTree(vector<BSPTree> &bsptree, int MAX_TREE_SIZE)
{
	for (int t = 0; t < bsptree.size(); t++)
	{
		bsptree[t].node = vector<vector<double>>(MAX_TREE_SIZE, vector<double>(DIMENSION));
		bsptree[t].fitness = vector<double>(MAX_TREE_SIZE, DBL_MAX);

		bsptree[t].axisline = vector<int>(MAX_TREE_SIZE, -1);
		bsptree[t].axisline[0] = -2;

		bsptree[t].min = vector<double>(MAX_TREE_SIZE, DBL_MAX);
		bsptree[t].max = vector<double>(MAX_TREE_SIZE, DBL_MAX);

		bsptree[t].sub_h = vector<vector<double>>(2, vector<double>(DIMENSION));
		bsptree[t].p_idx = vector<int>(MAX_TREE_SIZE, MAX_TREE_SIZE + 1);
		bsptree[t].l_idx = vector<int>(MAX_TREE_SIZE, MAX_TREE_SIZE + 1);
		bsptree[t].r_idx = vector<int>(MAX_TREE_SIZE, MAX_TREE_SIZE + 1);

		bsptree[t].cur_tree_size = 0;
	}
}

void SetInterval(vector<vector<double>> &sub_h, double lb, double ub)
{
	for (int j = 0; j < DIMENSION; ++j)
	{
		sub_h[0][j] = lb;
		sub_h[1][j] = ub;
	}
}

int SearchNode(BSPTree &cur_tree, vector<double> cur_data, double lb, double ub)
{
	if (cur_tree.axisline[0] == -2)
	{
		return -1;
	}
	else
	{
		int idx = 0;
		SetInterval(cur_tree.sub_h, lb, ub);
		while (cur_tree.axisline[idx] > -1)
		{
			int d = cur_tree.axisline[idx];
			if (cur_data[d] < (cur_tree.min[idx] + cur_tree.max[idx]) / 2.0)
			{
				cur_tree.sub_h[1][d] = cur_tree.max[idx];
				idx = cur_tree.l_idx[idx];
			}
			else
			{
				cur_tree.sub_h[0][d] = cur_tree.min[idx];
				idx = cur_tree.r_idx[idx];
			}
		}
		return idx;
	}
}

void InsertNode(BSPTree &cur_tree, int idx, vector<double> cur_data, double fitness)
{
	if (cur_tree.axisline[0] == -2) // ONLY ROOT WITH NOTHING
	{
		cur_tree.axisline[0] = -1;
		cur_tree.node[0] = cur_data;
		cur_tree.fitness[0] = fitness;
	}
	else
	{
		vector<double> memo_x = cur_tree.node[idx];
		double memo_f = cur_tree.fitness[idx];

		int d = 0;
		double d_diff = 0.0;
		for (int i = 0; i < cur_data.size(); i++)
		{
			double temp_diff = fabs(memo_x[i] - cur_data[i]);
			if (d_diff < temp_diff)
			{
				d = i;
				d_diff = temp_diff;
			}
		}

		cur_tree.axisline[idx] = d;
		cur_tree.min[idx] = min(memo_x[d], cur_data[d]);
		cur_tree.max[idx] = max(memo_x[d], cur_data[d]);

		int l_idx = cur_tree.cur_tree_size + 1;
		int r_idx = cur_tree.cur_tree_size + 2;
		cur_tree.cur_tree_size += 2;

		cur_tree.l_idx[idx] = l_idx;
		cur_tree.r_idx[idx] = r_idx;
		cur_tree.p_idx[l_idx] = idx;
		cur_tree.p_idx[r_idx] = idx;

		if (cur_data[d] < (cur_tree.min[idx] + cur_tree.max[idx]) / 2.0)
		{
			cur_tree.node[l_idx] = cur_data;
			cur_tree.fitness[l_idx] = fitness;
			cur_tree.node[r_idx] = memo_x;
			cur_tree.fitness[r_idx] = memo_f;
		}
		else
		{
			cur_tree.node[l_idx] = memo_x;
			cur_tree.fitness[l_idx] = memo_f;
			cur_tree.node[r_idx] = cur_data;
			cur_tree.fitness[r_idx] = fitness;
		}
	}
}

double CalculateDistance(vector<double> vec1, vector<double> vec2)
{
	double sigma = 0.0;
	for (int i = 0; i < vec1.size(); ++i)
		sigma += pow((vec1[i] - vec2[i]), 2);
	return sqrt(sigma);
}

void FreeNode(vector<BSPTree> &bsptree)
{
	bsptree.clear();
}
