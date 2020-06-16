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

int main(int argc, char* argv[])
{
	int cm_n=atoi(argv[1]);
	const int n = cm_n+1; //number of nodes + 1
	const int r = 1;

	double x[n];
	double y[n];

	x[0] = 5000; //200
	y[0] = 5000; //200

	for (int i = 1; i<n; i++)
	{
		x[i] = 0;
		y[i] = 0;
	}

	int node_generator_flag = 0;
	int min_pos_x = x[0]-r;
	int max_pos_x = x[0]+2*r;
	int min_pos_y = y[0]-r;
	int max_pos_y = y[0]+2*r;
	int seed=0;
	for (int i = 1; i < n; i++)
	{
	seed = 0;
	while (node_generator_flag == 0)
	{
	srand(time(NULL) + seed * 15);
	x[i] = rand() % (max_pos_x - min_pos_x) + min_pos_x;
	srand(time(NULL) + seed * 27);
	y[i] = rand() % (max_pos_y - min_pos_y) + min_pos_y;
	for (int j = 0; j < i; j++)
	if (distance(x[i], x[j], y[i], y[j]) == r)
	node_generator_flag = 1;
	for (int j = 0; j < i; j++)
	if (distance(x[i], x[j], y[i], y[j]) == 0)
	node_generator_flag = 0;
	seed++;
	}
	node_generator_flag = 0;
	if (x[i] >= max_pos_x)
	max_pos_x = x[i] + 2 * r;
	if (x[i] <= min_pos_x)
	min_pos_x = x[i] - r;
	if (y[i] >= max_pos_y)
	max_pos_y = y[i] + 2 * r;
	if (y[i] <= min_pos_y)
	min_pos_y = y[i] - r;
	}

	/*x[0] = 5000;
	y[0] = 5000;
	x[1] = 5000;
	y[1] = 5001;
	x[2] = 5001;
	y[2] = 5001;
	x[3] = 5000;
	y[3] = 4999;
	x[4] = 4999;
	y[4] = 5000;
	x[5] = 5000;
	y[5] = 4998;
	x[6] = 4999;
	y[6] = 4999;
	x[7] = 5001;
	y[7] = 4998;
	x[8] = 4999;
	y[8] = 4998;
	x[9] = 5001;
	y[9] = 4999;
	x[10] = 5001;
	y[10] = 4997;
	x[11] = 4998;
	y[11] = 4999;
	x[12] = 5001;
	y[12] = 4996;
	x[13] = 4998;
	y[13] = 5000;
	x[14] = 4998;
	y[14] = 5001;*/

	// preparing results for colormap picture /////////////////////////////////////////////////////////////////////////////////////////////

	// writing x to file
	ofstream file_x;
	file_x.open("S04_x.txt");
	for(int i=1; i<n; i++)
		file_x<<i<<" "<<x[i]<<"\n";
	file_x<<n<<" "<<x[0]<<"\n";
	file_x.close();

	// writing y to file
	ofstream file_y;
	file_y.open("S05_y.txt");
	for(int i=1; i<n; i++)
		file_y<<i<<" "<<y[i]<<"\n";
	file_y<<n<<" "<<y[0]<<"\n";
	file_y.close();

	// writing set N to file
	ofstream file_set_N;
	file_set_N.open("S01_set_N.txt");
	for (int i = 1; i < n; i++)
	{
		file_set_N << i << "\n";
	}
	file_set_N.close();

	// writing set A to file
	ofstream file_set_A;
	file_set_A.open("S01_set_A.txt");
	for (int i = 1; i <= n; i++)
	{
		file_set_A << i << "\n";
	}
	file_set_A.close();

	// writing Location of Nodes to file
	ofstream file_LoN;
	file_LoN.open("R03_LoN.txt", std::ios::app);
	for (int i = 0; i < n; i++)
	{
		file_LoN<<"x["<<i<<"] = "<<x[i]<<";"<<"\n";
		file_LoN<<"y["<<i<<"] = "<<y[i]<<";"<<"\n";
	}
	file_LoN.close();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}
