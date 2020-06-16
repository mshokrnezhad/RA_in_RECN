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

			for (int k2 = 1; k2 < no_specs; k2++)
			{
				for (int i = 0; i < n - 1; i++)
					for (int j = 0; j < n; j++)
						 a[i][j] = 0;

				for (int i = 0; i < n - 1; i++)
					for (int j = 0; j < n; j++)
						if (channel[i + 1] != k2)
						{
							if (i == j)
								a[i][j] = 1;
							else
								a[i][j] = 0;
						}
						else
						{
							if (channel[j + 1] != k2)
								a[i][j] = 0;
							if (channel[j + 1] == k2 && i == j)
								a[i][j] = 1;
							if (channel[j + 1] == k2 && i != j)
								a[i][j] = ((-target_sinr)*h(x[j + 1], x[next_node[i + 1]], y[j + 1], y[next_node[i + 1]])) / h(x[i + 1], x[next_node[i + 1]], y[i + 1], y[next_node[i + 1]]);
							a[i][n - 1] = ((target_sinr)*noise) / h(x[i + 1], x[next_node[i + 1]], y[i + 1], y[next_node[i + 1]]);
						}

				int i, j, k;
				float xge[n - 1];
				int size = n - 1;

				for (i = 0; i < size; i++)                    //Pivotisation
					for (k = i + 1; k < size; k++)
						if (a[i][i] < a[k][i])
							for (j = 0; j <= size; j++)
							{
								double tge = a[i][j];
								a[i][j] = a[k][j];
								a[k][j] = tge;
							}

				for (i = 0; i < size - 1; i++)            //loop to perform the gauss elimination
					for (k = i + 1; k < size; k++)
					{
						double t = a[k][i] / a[i][i];
						for (j = 0; j <= size; j++)
							a[k][j] = a[k][j] - t*a[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
					}

				for (i = size - 1; i >= 0; i--)                //back-substitution
				{                        //x is an array whose values correspond to the values of x,y,z..
					xge[i] = a[i][size];                //make the variable to be calculated equal to the rhs of the last equation
					for (j = 0; j < size; j++)
						if (j != i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
							xge[i] = xge[i] - a[i][j] * xge[j];
					xge[i] = xge[i] / a[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
				}

				for (int i = 0; i < size; i++)
					if (xge[i] < 0)
						xge[i] = max_power;
				double tge = 0;

				for (i = 0; i < size; i++)
					tge = xge[i] + tge;            // Print the values of x, y,z,....    

				tsop[k1] = tge + tsop[k1];
			}
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

	// writing execution time to file
	ofstream file_et;
	file_et.open("R04_cpp_chat.txt", std::ios::app);
	file_et << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
	file_et.close();
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}
