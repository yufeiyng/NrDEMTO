#pragma once
#include <vector>

using namespace std;

struct BSPTree
{
	vector<vector<double>> node;
	vector<double> fitness;

	vector<int> axisline;
	vector<vector<double>> sub_h;
	vector<double> min;
	vector<double> max;

	vector<int> p_idx;
	vector<int> l_idx;
	vector<int> r_idx;

	int cur_tree_size;
};

void InitializeTree(vector<BSPTree> &bsptree, int MAX_TREE_SIZE);

void SetInterval(vector<vector<double>> &sub_h, double lb, double ub);

int SearchNode(BSPTree &cur_tree, vector<double> cur_data, double lb, double ub);

void InsertNode(BSPTree &cur_tree, int idx, vector<double> cur_data, double fitness);

double CalculateDistance(vector<double> vec1, vector<double> vec2);

void FreeNode(vector<BSPTree> &bsptree);