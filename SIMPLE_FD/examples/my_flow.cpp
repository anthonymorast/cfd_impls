#include <iostream>

using namespace std;

#define EPS 0.0000001     // convergence for pressure
#define COVMAX 0.000015   // psi iteration
#define KMAX 4200         // ???
#define ALPHA 0.6f        // underrelaxation factor


double** alloc_table(int nx, int ny);
void apply_boundaries(double**, double**, double**) ;


int main(int argc, char* argv[])
{
  int nx = 50;
  int ny = 50;

  double dx = (1.0f/(double)nx);
  double dy = (1.0f/(double)ny);
  double dt = 0.01;

  double lid_vel = 1.0f;
  double re = 100;

  double **u = alloc_table(nx-1, ny);  // [NX+1][NY]
  double **v = alloc_table(nx, ny-1);  // [NX][NY+1]
  double **p_hat = alloc_table(nx, ny);
  double **p_star = alloc_table(nx, ny);
  double **u_star = alloc_table(nx-1, ny);
  double **v_star = alloc_table(nx, ny-1);
  double **p = alloc_table(nx, ny);  // [NX][NY]

  double cov;


}

void apply_boundaries(double **p, double **u, double **v, int nx, int ny)
{
  // pressure boudnaries -- dp/dt = 0
  for(int i = 0; i < nx; i++)
  {
    p[i][0] = p[i][1];        // lid
    p[i][ny-1] = p[i][ny-2];  // cavity bottom
  }

/*  for (int i = 0; i < ny; i++)
  {
    p[0][j] = p[1][j];        // left wall (inlet for channel)
    p[nx-1][j] = p[nx-2][j]   // right wall (outlet for channel)
  } */


}

double** alloc_table(int nx, int ny)
{
  double **tab = new double*[nx];
  for(int i = 0; i < nx; i++)
  {
      tab[i] = new double[ny];
  }

  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < ny; j++)
    {
      tab[i][j] = 0;
    }
  }

  return tab;
}
