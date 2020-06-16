// With noncooperative power update function
//
#include "iostream"
#include "fstream"
#include "iomanip"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "vector"
#include "ctime"

using namespace std;

double distance(double x_i, double x_j, double y_i, double y_j)
{
	return sqrt(pow((x_i - x_j), 2) + pow((y_i - y_j), 2));
}

double h(double x_i, double x_next_node_j, double y_i, double y_next_node_j)
{
	if (distance(x_i, x_next_node_j, y_i, y_next_node_j) == 0)
		return 1;
	else
		return 0.09*pow(distance(x_i, x_next_node_j, y_i, y_next_node_j), -3);
}

int random_generator(int min, int max)
{
	int random_number;
	//srand(time(NULL));
	random_number = rand() % (1000 - 0) + 0;
	for (int i = 0; i<max - min + 1; i++)
		if (random_number >= ((1000 * i) / (max - min + 1)) && random_number <= ((1000 * (i + 1)) / (max - min + 1)))
			return i + min;
}

void print_progress_bar(int percent)
{
	string bar;
	cout << "\r" << bar;
	cout << percent << "% " << std::flush;
}

int func_next_node(int i, int n, int r, double x[], double y[])
{
	int temp_distance = 1000000;
	int n_n = 0;
	if (i == 0)
		return n_n;
	else
	{
		for (int j = 0; j < i; j++)
			if (i != j)
				if (distance(x[i], x[j], y[i], y[j]) <= r)
					if (distance(x[j], x[0], y[j], y[0]) <= temp_distance)
					{
						temp_distance = distance(x[j], x[0], y[j], y[0]);
						n_n = j;
					}
		return n_n;
	}
}

double func_get_max_index(double arr[], int size)
{
	int MaxIndex;
	double temp_max = 0;
	for (int i = 0; i<size; i++)
		if (arr[i]>temp_max)
		{
			temp_max = arr[i];
			MaxIndex = i;
		}

	return MaxIndex;
}

double func_get_min_index(double arr[], int size)
{
	int MinIndex;
	double temp_min = 1000000000000;
	for (int i = 0; i<size; i++)
		if (arr[i]<temp_min)
		{
			temp_min = arr[i];
			MinIndex = i;
		}

	return MinIndex;
}

double func_get_max(double arr[], int size)
{
	int MaxIndex;
	double temp_max = 0;
	for (int i = 0; i<size; i++)
		if (arr[i]>temp_max)
		{
			temp_max = arr[i];
			MaxIndex = i;
		}

	return temp_max;
}

double func_get_min(double arr[], int size)
{
	int MinIndex;
	double temp_min = 1000000000000;
	for (int i = 0; i<size; i++)
		if (arr[i]<temp_min)
		{
			temp_min = arr[i];
			MinIndex = i;
		}

	return temp_min;
}

double diff(double a, double b)
{
	if (a >= b)
		return a - b;
	else
		return b - a;
}

int main(int argc, char* argv[])
{
	//int cm_n = atoi(argv[1]);
	//int cm_no_specs = atoi(argv[2]);
	//double cm_target_sinr = atof(argv[3]);
	const int n = 6 + 1;
	const int no_specs = 2 + 1;
	double target_sinr = 0.5; //0.05; //4;
	const int hop = 1;
	const int r = 1;
	double noise = 0.0000000001; //0.0000000001;
	double max_power = 2;
	const int round_bound = 20; //minimun amount is 1 where we have just one round
	const int iteration_bound = 500;

	double sinr[round_bound][n][no_specs];
	double t_sinr[iteration_bound][n][no_specs];
	double sum_sinr[round_bound];
	int sinr_received_counter = 0;

	double p[round_bound][n][no_specs];
	double t_p[iteration_bound][n][no_specs];
	double sum_power[round_bound];

	double I[round_bound][n][no_specs];
	double t_I[iteration_bound][n][no_specs];

	double outage_ratio[round_bound];

	int max_powered_node;
	int min_interfered_channel;

	int round_loop_index_counter = 0;
	int round_loop_flag = 0;

	int next_node[n];
	int chflag[n];
	int round;
	int channel[round_bound][n];
	int step_channel[n][n];
	double temp[n + no_specs];

	int iteration_loop_index_counter = 0;
	int iteration_loop_flag = 0;
	int iteration;

	double x[n];
	double y[n];

	double noh[n];
	int pgsan[n];
	int step[n];
	int s = 1;

	int temp_a = 0;
	int temp_c = 0;
	double temp_b = 0;

	ofstream file_dP;
	file_dP.open("R01_dgame_dP.txt", std::ios::app);
	ofstream file_dSIR;
	file_dSIR.open("R01_dgame_dSIR.txt", std::ios::app);

	ifstream rxfile;
	rxfile.open("S04_x.txt");
	while (!rxfile.eof())
	{
		rxfile >> temp_a >> temp_b;
		if (temp_a != n)
			x[temp_a] = temp_b;
		else
			x[0] = temp_b;
	}
	rxfile.close();

	temp_a = 0;
	temp_b = 0;
	ifstream ryfile;
	ryfile.open("S05_y.txt");
	while (!ryfile.eof())
	{
		ryfile >> temp_a >> temp_b;
		if (temp_a != n)
			y[temp_a] = temp_b;
		else
			y[0] = temp_b;
	}
	ryfile.close();

	int start_s = clock();

	/*for (int i = 0; i<n; i++)
	next_node[i] = func_next_node(i, n, r, x, y);*/
	//*****************************************************
	next_node[0] = 0;
	next_node[1] = 0;
	next_node[2] = 0;
	next_node[3] = 1;
	next_node[4] = 1;
	next_node[5] = 2;
	next_node[6] = 2;
	//*****************************************************

	for (int i = 0; i < n; i++)
		noh[i] = 0;

	for (int j = 1; j < n; j++)
		noh[j] = noh[next_node[j]] + 1;
	noh[0] = n;

	for (int i = 1; i < n; i++)
	{
		pgsan[i] = func_get_min_index(noh, n);
		noh[pgsan[i]] = n;
	}
	pgsan[0] = 0;

	step[0] = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 1; j < n; j++)
			if (next_node[j] == pgsan[i])
				step[j] = s;
		s++;
	}

	for (int i = 1; i < n; i++)
		for (int r = 0; r < round_bound; r++)
		{
			channel[r][i] = 0;// temp_nu%no_specs;
			chflag[i] = 0;
		}

	for (int s = 1; s<n; s++)
	{
		round = 0;
		round_loop_flag = 0;
		//while (sinr_received_counter != no_trans && round<round_bound && round_loop_flag != 1)
		while (round < round_bound && round_loop_flag != 1)
		{
			file_dP << "s" << s << "r" << round << "\n";
			file_dSIR << "s" << s << "r" << round << "\n";

			// first assignment of channel and power to all nodes
			if (round == 0)
				for (int i = 1; i < n; i++)
				{
					if (chflag[i] == 1)
						channel[round][i] = step_channel[s - 1][i];// temp_nu%no_specs;
					if (step[i] == s)
						channel[round][i] = 1;// temp_nu%no_specs;
				}
			else
			{
				// finding max_interfered_node

				for (int i = 0; i<n + no_specs; i++)
					temp[i] = 0;

				for (int i = 1; i < n; i++)
					if (chflag[i] == 0)
						temp[i] = p[round - 1][i][channel[round - 1][i]];
					else
						temp[i] = 0;
				max_powered_node = func_get_max_index(temp, n);

				// finding min_interfered_channel
				for (int k = 1; k<no_specs; k++)
					temp[k - 1] = I[round - 1][max_powered_node][k];
				min_interfered_channel = func_get_min_index(temp, no_specs - 1) + 1;

				for (int i = 1; i<n; i++)
				{
					// assigning channel and power to max_interfered_node
					if (i == max_powered_node)
						channel[round][i] = min_interfered_channel;

					// assigning channel and power to other nodes
					else
						channel[round][i] = channel[round - 1][i];
				}
			}

			// power control procedure
			iteration = 0;
			for (int i = 1; i < n; i++)
			{
				for (int k = 0; k < no_specs; k++)
					t_p[iteration][i][k] = 0;
				if(channel[round][i]==0)
					t_p[iteration][i][channel[round][i]] = 0;
				else
					t_p[iteration][i][channel[round][i]] = max_power;
			}

			iteration_loop_flag = 0;

			// power control itterations
			while (iteration < iteration_bound && iteration_loop_flag != 1)
			{
				if (iteration != 0)
					for (int i = 1; i < n; i++)
					{
						for (int k = 1; k < no_specs; k++)
							t_p[iteration][i][k] = 0;
						//p[iteration][i][channel[round][i]] = target_sinr*(p[iteration - 1][i][channel[round][i]] / sinr[iteration - 1][i][channel[round][i]]);
						if(channel[round][i]==0)
							t_p[iteration][i][channel[round][i]] = 0;
						else
							t_p[iteration][i][channel[round][i]] = min(max_power, target_sinr*(t_p[iteration - 1][i][channel[round][i]] / t_sinr[iteration - 1][i][channel[round][i]]));
					}

				for (int i = 1; i < n; i++)
					for (int k = 0; k < no_specs; k++)
					{
						t_I[iteration][i][k] = 0;
						for (int j = 1; j < n; j++)
							if (j != i)
								t_I[iteration][i][k] = t_I[iteration][i][k] + t_p[iteration][j][k] * h(x[j], x[next_node[i]], y[j], y[next_node[i]]);
						t_I[iteration][i][k] = t_I[iteration][i][k] + noise;
						t_sinr[iteration][i][k] = (t_p[iteration][i][k] * h(x[i], x[next_node[i]], y[i], y[next_node[i]])) / t_I[iteration][i][k];
					}
				for (int i = 1; i < n; i++)
				{
					file_dP << i << " " << t_p[iteration][i][channel[round][i]] << "\n";
					file_dSIR << i << " " << t_sinr[iteration][i][channel[round][i]] << "\n";
				}

				iteration++;

				iteration_loop_index_counter = 0;

				if (iteration > 5)
				{
					for (int it = iteration - 2; it >= iteration - 5; it--)
						for (int i = 1; i < n; i++)
							if (diff(t_sinr[iteration - 1][i][channel[round][i]], t_sinr[it][i][channel[round][i]]) == 0)
								iteration_loop_index_counter++;
					if (iteration_loop_index_counter == 4 * (n - 1))
						iteration_loop_flag = 1;
				}
			}

			for (int i = 1; i < n; i++)
				for (int k = 0; k < no_specs; k++)
				{
					p[round][i][k] = t_p[iteration - 1][i][k];
					I[round][i][k] = t_I[iteration - 1][i][k];
					sinr[round][i][k] = t_sinr[iteration - 1][i][k];
				}
			
			//**********************************************************************************************************************************
			//Without this line, the program was executed in MS. VS	perfectly but in Linux this part should be added to produce flawless results		
			for (int i = 1; i < n; i++)
				p[round][i][0] = 0;
			//**********************************************************************************************************************************

				// final values of each power control procedure
				sinr_received_counter = 0;
				sum_sinr[round] = 0;
				sum_power[round] = 0;
				for (int i = 1; i<n; i++)
				{
					if (diff(sinr[round][i][channel[round][i]], target_sinr)<0.01)
						sinr_received_counter++;

					for (int k = 1; k<no_specs; k++)
					{
						sum_sinr[round] = sum_sinr[round] + sinr[round][i][k];
						sum_power[round] = sum_power[round] + p[round][i][k];
					}
				}
				outage_ratio[round] = 1 - (sinr_received_counter / (double)(n - 1));

				if (round > 0)
				{
					if (sum_power[round] > sum_power[round - 1])
						for (int i = 1; i < n; i++)
						{
							channel[round][i] = channel[round - 1][i];
							for (int k = 0; k < no_specs; k++)
							{
								I[round][i][k] = I[round - 1][i][k];
								p[round][i][k] = p[round - 1][i][k];
								sum_power[round] = sum_power[round - 1];
								sinr[round][i][k] = sinr[round - 1][i][k];
							}
						}
					else
						chflag[max_powered_node] = 1;
				}

				round_loop_flag = 0;
				round_loop_index_counter = 0;
				for (int i = 1; i < n; i++)
					if (step[i] == s && chflag[i] == 1)
						round_loop_index_counter++;
				for (int i = 1; i < n; i++)
					if (step[i] == s)
						round_loop_index_counter--;
				if (round_loop_index_counter == 0)
					round_loop_flag = 1;


				/*cout << "Channel:   ";
				for (int i = 1; i<n; i++)
				cout << channel[round][i] << " ";
				cout << "\n";
				cout << "Path Gain: ";
				for (int i = 1; i<n; i++)
				cout << h(x[i], x[next_node[i]], y[i], y[next_node[i]]) << " ";
				cout << "\n";
				cout << "Power:     ";
				for (int i = 1; i<n; i++)
				cout << p[round][i][channel[round][i]] << " ";
				cout << "\n";
				cout << "SINR:      ";
				for (int i = 1; i<n; i++)
				cout << sinr[round][i][channel[round][i]] << " ";
				cout << "\n";
				cout << sum_power[round] << " " << sum_sinr[round] << "\n";*/

				file_dP << "s" << s << "r" << round << "\n";
				file_dSIR << "s" << s << "r" << round << "\n";


				round++;
			}
			for (int i = 1; i < n; i++)
				step_channel[s][i] = channel[round - 1][i];
		}

		int stop_s = clock();

		file_dP.close();
		file_dSIR.close();

		// preparing results for colormap picture /////////////////////////////////////////////////////////////////////////////////////////////

		/*

		// writing Sum of Power to file
		ofstream file_SoP;
		file_SoP.open("R01_game_SoP.txt", std::ios::app);
		file_SoP << sum_power[round - 1] << "\n";
		file_SoP.close();

		// writing Sum of SINR to file
		ofstream file_SoS;
		file_SoS.open("R02_game_SoS.txt", std::ios::app);
		file_SoS << sum_sinr[round - 1] << "\n";
		file_SoS.close();

		// writing execution time to file
		ofstream file_et;
		file_et.open("R04_game_t.txt", std::ios::app);
		file_et << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
		file_et.close();

		*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}
