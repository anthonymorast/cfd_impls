#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

void dynamic_memcpy(double **a, double **b, int rows, int cols);
double** alloc_table(int numx, int numy);   // initialzies the variables to 0
void print_array(double **arr, int x_max, int y_max);
void array_to_csv(double **arr, int numx, int numy, string filename, string var, int iter);
