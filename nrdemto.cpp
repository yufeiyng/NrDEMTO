#include "de.h"
#include "bsptree.h"

#include <stdlib.h>
#include <time.h>
#include <random>
#include <fstream>
#include <iostream>

int main()
{
	srand((unsigned int)time(NULL));

	std::default_random_engine engine;
	std::cauchy_distribution<double> f(0.3, 0.1);
	std::normal_distribution<double> cr(0.4, 0.1);

	std::ofstream ofs_t1;
	std::ofstream ofs_t2;
	ofs_t1.open("nrdemto_t1.csv", std::ios::out);
	ofs_t2.open("nrdemto_t2.csv", std::ios::out);

	int num_of_task;
	int MAX_TREE_SIZE = 2 * MAX_EVALUATIONS;
	
	std::cout << "PLEASE ENTER THE NUMBER OF TASK(1~9):\n";
	std::cin >> num_of_task;

	clock_t start = clock();

	for (int r = 0; r < MAX_RUN_TIMES; r++)
	{
		vector<int> cur_evals(T, 0);
		vector<BSPTree> bsp(T);
		vector<DE> de(T);

		InitializeTree(bsp, MAX_TREE_SIZE);
		InitializeDE(de);
		SetBound(de, num_of_task);
		for (int t = 0; t < T; ++t)
		{
			InitializePop(de[t]);
			EvaluatePop(bsp[t], de, de[t].pop, de[t].fitness, cur_evals, t, num_of_task);
		}

		while (cur_evals[0] + cur_evals[1] < MAX_EVALUATIONS)
		{
			for (int t = 0; t < T; t++)
			{
				Evolution(de[t], f(engine), cr(engine), t);
				EvaluatePop(bsp[t], de, de[t].trial, de[t].f_trial, cur_evals, t, num_of_task);
				Selection(de[t]);
				// TRANSFER PROCESS
				if (GetRandom() < RMP)
				{
					Transfer(de, t);
				}

				// LEARNING SCHEME BASED ON CURRENT BEST SOLUTION
				int idx = SearchNode(bsp[t], de[t].best_solution, de[t].lb, de[t].ub);
				int i_rand = rand() % POP_SIZE;
				LocalSearch(bsp[t], idx, de, i_rand, t, cur_evals, num_of_task);

				/*if (c % 10 == 0)
				{
					if (t == 0)
						ofs_t1 << de[t].best_fitness << std::endl;
					else
						ofs_t2 << de[t].best_fitness << std::endl;
				}*/
			}
		}

		for (int t = 0; t < T; t++)
		{
			printf("%.4lE\n", de[t].best_fitness);
			// BOX_PLOT
			if (t == 0)
				ofs_t1 << de[t].best_fitness << std::endl;
			else
				ofs_t2 << de[t].best_fitness << std::endl;
		}
		std::cout << std::endl;

		bsp.clear();
		de.clear();
		cur_evals.clear();
	}

	clock_t stop = clock();
	float cpu_time = (float)(stop - start);
	printf("nrdemto time (CPU): %f ms \n", cpu_time);

	ofs_t1.close();
	ofs_t2.close();

	return 0;
}