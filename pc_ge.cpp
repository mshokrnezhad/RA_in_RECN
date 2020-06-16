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
	const int n = cm_n+1;
	const int no_specs = cm_no_specs+1;
	double target_sinr = cm_target_sinr; //0.05; //4;
	const int hop = 1;
	const int r = 1;
	double noise = 1; //0.0000000001;
	double max_power = 10000;
	


	double p[n];
	double a[n - 1][n];
	
	double sinr[n];
	double sum_sinr;
	int sinr_received_counter = 0;

	double sum_power;

	double I[n];

	double outage_ratio;

	int next_node[n];
	int channel[n];

	double x[n];
	double y[n];

	for (int i = 0; i<n; i++)
		channel[i] = 0;

	int temp_a = 0;
	int temp_c = 0;
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

	temp_a = 0;
	temp_b = 0;
	ifstream chfile;
	chfile.open("S06_cpp_ch.txt");
	while (!chfile.eof())
	{
		chfile >> temp_a >> temp_c >> temp_b;
		if (temp_b > 0.9)
			channel[temp_a] = temp_c;
	}
	chfile.close();

	temp_a = 0;
	temp_b = 0;
	ifstream nnfile;
	nnfile.open("S03_nn.txt");
	while (!nnfile.eof())
	{
		nnfile >> temp_a >> temp_b;
		if (temp_b == n)
			next_node[temp_a] = 0;
		else
			next_node[temp_a] = temp_b;
	}
	nnfile.close();
	next_node[0]=0;

	int start_s = clock();


	for (int i = 0; i < n; i++)
		p[i] = 0;

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
			if (xge[i] < 0 )
				xge[i] = max_power;

		for (int i = 0; i < size; i++)
			if (xge[i]>0)
				p[i + 1] = xge[i];
	}

	for (int i = 1; i < n; i++)
	{
		I[i] = 0;
		for (int j = 1; j < n; j++)
			if (j != i && channel[i] == channel[j])
				I[i] = I[i] + p[j] * h(x[j], x[next_node[i]], y[j], y[next_node[i]]);
		I[i] = I[i] + noise;
		sinr[i] = (p[i] * h(x[i], x[next_node[i]], y[i], y[next_node[i]])) / I[i];
	}

	// final values of each power control procedure
	sinr_received_counter = 0;
	sum_sinr = 0;
	sum_power = 0;
	for (int i = 1; i<n; i++)
	{
		if (sinr[i] == target_sinr)
			sinr_received_counter++;

		sum_sinr = sum_sinr + sinr[i];
		sum_power = sum_power + p[i];
	}
	outage_ratio = 1 - (sinr_received_counter / (double)(n - 1));

	/*cout << "Channel:   ";
	for (int i = 1; i<n; i++)
		cout << channel[i] << " ";
	cout << "\n";
	cout << "Path Gain: ";
	for (int i = 1; i<n; i++)
		cout << h(x[i], x[next_node[i]], y[i], y[next_node[i]]) << " ";
	cout << "\n";
	cout << "Power:     ";
	for (int i = 1; i<n; i++)
		cout << p[i] << " ";
	cout << "\n";
	cout << "SINR:      ";
	for (int i = 1; i<n; i++)
		cout << sinr[i] << " ";
	cout << "\n";
	cout << sum_power << " " << sum_sinr << "\n";*/

	int stop_s = clock();


	// writing Sum of Power to file
	ofstream file_SoP;
	file_SoP.open("R01_cpp_SoP.txt", std::ios::app);
	file_SoP << sum_power << "\n";
	file_SoP.close();

	// writing Sum of SINR to file
	ofstream file_SoS;
	file_SoS.open("R02_cpp_SoS.txt", std::ios::app);
	file_SoS << sum_sinr << "\n";
	file_SoS.close();

	// writing execution time to file
	ofstream file_et;
	file_et.open("R04_cpp_pct.txt", std::ios::app);
	file_et << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
	file_et.close();

	// reporting results 

	/*cout << fixed << setprecision(15);
	cout << "#" << no_specs << "_" << n << "_" << 1 << " " << max_sum_sinr << " " << min_sum_power << " " << min_outage_ratio << " " << no_trans << "\n";*/


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}
