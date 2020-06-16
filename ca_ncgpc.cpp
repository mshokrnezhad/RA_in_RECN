// iot_pc_cha_d2a.cpp : Defines the entry point for the console application.
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
	int cm_n=atoi(argv[1]);
	int cm_no_specs=atoi(argv[2]);
	double cm_target_sinr=atof(argv[3]);

	const int n = cm_n + 1;
	const int no_specs = cm_no_specs + 1;
	double target_sinr = cm_target_sinr; //0.05; //4;
	const int hop = 1;
	const int r = 1;
	double noise = 1; //0.0000000001;
	double max_power = 10000;
	//double opc_constant = 0.01;

	int next_node[n];
	//double temp[n + no_specs];
	int pgsan[n - 1]; //min to max path-gain sorted arrey of nodes
	int channel[n];
	//int nuch[no_specs]; //number of nodes in each channel
	//double spgch[no_specs]; //sume of path-gains in each channel
	double noh[n]; //devide i's number of hops to BS
	double a[n - 1][n];
	double tsop[no_specs];

	const int iteration_bound = 200;

	double t_sinr[iteration_bound][n][no_specs];
	double t_p[iteration_bound][n][no_specs];
	double t_I[iteration_bound][n][no_specs];
	int iteration = 0;

	for (int i = 0; i < n; i++)
	{
		if (i < n - 1)
			pgsan[i] = 0;
		channel[i] = 0;
		noh[i] = 0;
	}

	double x[n];
	double y[n];
	int temp_a = 0;
	double temp_b = 0;

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

	for (int i = 0; i < n; i++)
		next_node[i] = func_next_node(i, n, r, x, y);

	for (int j = 1; j < n; j++)
		noh[j] = noh[next_node[j]] + 1;
	noh[0] = n;

	for (int i = 0; i < n - 1; i++)
	{
		pgsan[i] = func_get_min_index(noh, n);
		noh[pgsan[i]] = n;
	}

	for (int i = 1; i < n; i++)
		channel[i] = 0;

	for (int round = 0; round < n - 1; round++)
	{
		for (int k1 = 1; k1 < no_specs; k1++)
		{
			channel[pgsan[round]] = k1;
			tsop[k1] = 0;

			// power control procedure
			iteration = 0;
			for (int i = 1; i < n; i++)
			{
				for (int k2 = 0; k2 < no_specs; k2++)
					t_p[iteration][i][k2] = 0;
				if(channel[i]==0)
					t_p[iteration][i][channel[i]] = 0;
				else
					t_p[iteration][i][channel[i]] = max_power;
			}

			// power control itterations
			while (iteration < iteration_bound)
			{
				if (iteration != 0)
					for (int i = 1; i < n; i++)
					{
						for (int k3 = 1; k3 < no_specs; k3++)
							t_p[iteration][i][k3] = 0;
						//p[iteration][i][channel[round][i]] = target_sinr*(p[iteration - 1][i][channel[round][i]] / sinr[iteration - 1][i][channel[round][i]]);
						if(channel[i]==0)
							t_p[iteration][i][channel[i]] = 0;
						else
							t_p[iteration][i][channel[i]] = min(max_power, target_sinr*(t_p[iteration - 1][i][channel[i]] / t_sinr[iteration - 1][i][channel[i]]));
					}

				for (int i = 1; i < n; i++)
					for (int k4 = 0; k4 < no_specs; k4++)
					{
						t_I[iteration][i][k4] = 0;
						for (int j = 1; j < n; j++)
							if (j != i)
								t_I[iteration][i][k4] = t_I[iteration][i][k4] + t_p[iteration][j][k4] * h(x[j], x[next_node[i]], y[j], y[next_node[i]]);
						t_I[iteration][i][k4] = t_I[iteration][i][k4] + noise;
						t_sinr[iteration][i][k4] = (t_p[iteration][i][k4] * h(x[i], x[next_node[i]], y[i], y[next_node[i]])) / t_I[iteration][i][k4];
					}

				iteration++;
			}


			for (int i = 1; i < n; i++)
				for (int k = 0; k < no_specs; k++)
					tsop[k1]=tsop[k1]+t_p[iteration - 1][i][k];

		}

		tsop[0] = func_get_max(tsop, no_specs) + 1000;
		channel[pgsan[round]] = func_get_min_index(tsop, no_specs);
	}

	int stop_s = clock();	
	
	// writing assigned channels to the file
	ofstream file_channels;
	file_channels.open("S06_cpp_ch.txt");
	for (int i = 1; i < n; i++)
	file_channels << i << " " << channel[i] << " "<< "1" << "\n";
	file_channels.close();

	// writing next_node array to file
	ofstream file_next_node;
	file_next_node.open("S03_nn.txt");

	for (int i = 1; i < n; i++)
	{
	if (next_node[i] == 0)
	{
	file_next_node << i << " " << n << "\n";
	}
	else
	file_next_node << i << " " << next_node[i] << "\n";
	}
	file_next_node.close();

	/*// writing set K to file
	ofstream file_set_K;
	file_set_K.open("S02_set_K.txt");
	for (int i = 1; i < no_specs; i++)
	{
	file_set_K << i << "\n";
	}
	file_set_K.close();*/

	/*// writing execution time to file
	ofstream file_et;
	file_et.open("R04_cpp_chat.txt", std::ios::app);
	file_et << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
	file_et.close();*/
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}
