#include "utils.h"

/*
*   Copies array a to array b.
*/
void dynamic_memcpy(double **a, double **b, int rows, int cols)
{
    for(int i = 0; i < rows; i++)
    {
      memcpy(b[i], a[i], cols*sizeof(double));
    }
}

/*
*   Allocates a 2D array of doubles, initializing the array to 0's
*/
double** alloc_table(int numx, int numy)
{
    double **tab = new double*[numy];
    for(int i = 0; i < numy; i++)
    {
        tab[i] = new double[numx];
        for(int j = 0; j < numx; j++)
            tab[i][j] = 0;
    }

    return tab;
}

/*
*   Prints a 2D array of doubles as though it were a grid
*/
void print_array(double **arr, int x_max, int y_max)
{
    for(int i = 0; i < y_max; i++)
    {
        for(int j = 0; j < x_max; j++)
        {
            cout << arr[i][j] << ", ";
        }
        cout << endl;
    }
}

/*
*   Output a 2D array to a CSV file. (Could probably combine this with print_array()
*   by passing in the output stream (COUT or OFSTREAM)).
*/
void array_to_csv(double **arr, int numx, int numy, string filename, string var, int iter)
{
    ofstream out;
    out.open(filename, fstream::out);
    if(!out.is_open())
    {
        cout << "Error opening file '" << filename << "'." << endl;
        return;
    }

    out << "x,y,t," << var << endl;
    for(int i = 0; i < numy; i++)
    {
        for(int j = 0; j < numx; j++)
        {
            out << i << "," << j << "," << iter << "," << arr[i][j] << endl;
        }
    }
}
