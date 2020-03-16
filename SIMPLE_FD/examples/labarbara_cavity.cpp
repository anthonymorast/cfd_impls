#include <iostream>
#include <cstring>

using namespace std;

/* globals */
int NX = 51;
int NY = 51;

/* function defs */
void pressure_rhs(double **b, double **u, double **v, double rho, double dt,
  double dx, double dy);
void pressure_poisson(double **p, double dx, double dy, double** b, int max_iter, double eps);
void dynamic_memcpy(double **a, double **b, int rows, int cols);
double** alloc_table(int numx, int numy);   // initialzies the variables to 0
void cavity_flow(double **u, double **v, double **p, double **p_rhs, double density,
                 double viscosity, int nt, double dx, double dy, double dt,
                 int max_iter, double eps);

int main(int argc, char* argv[])
{
  // TODO: create a struct/class for flow params

  /* Poisson Pressure iteration variables */
  int max_iter = 100;    // max iterations for Poisson equation
  double epsilon = 0.0000001;

  double u_vel = 1.0;      // velocity at lid
  double max_x = 2;
  double max_y = 2;
  double dx = max_x / (NX - 1);
  double dy = max_y / (NY - 1);
  int nt = 500;

  // Re = lid length * lid velocity * density / viscosity
  //    = max_x * u_vel * density / viscosity
  //    = 2 * 1 * 1 / 0.005 =  400
  double density = 1.0;     // rho
  double viscosity = 0.005;   // nu
  double dt = 0.001;

  double **u = alloc_table(NX, NY);
  double **v = alloc_table(NX, NY);
  double **p = alloc_table(NX, NY);
  double **p_rhs = alloc_table(NX, NY);

  cavity_flow(u, v, p, p_rhs, density, viscosity, nt, dx, dy, dt, max_iter, epsilon);

  for(int i = 0; i < NX; i++)
  {
    for(int j = 0; j < NY; j++)
    {
      cout << u[i][j] << " ";
    }
    cout << endl;
  }

  return 0;
}

void cavity_flow(double **u, double **v, double **p, double **p_rhs, double density,
                 double viscosity, int nt, double dx, double dy, double dt,
                 int max_iter, double eps)
{
  double** un = alloc_table(NX, NY);
  double** vn = alloc_table(NX, NY);

  for(int t = 0; t < nt; t++)
  {
    dynamic_memcpy(u, un, NX, NY);
    dynamic_memcpy(v, vn, NX, NY);

    pressure_rhs(p_rhs, u, v, density, dt, dx, dy);
    pressure_poisson(p, dx, dy, p_rhs, max_iter, eps);

    for(int i = 1; i < (NX-1); i++)
    {
      for(int j = 1; j < (NY-1); j++)
      {
        u[i][j] = (
            (un[i][j]) -
            (un[i][j] * (dt/dx) * (un[i][j] - un[i-1][j])) -
            (vn[i][j] * (dt/dy) * (un[i][j] - un[i][j-1])) -
            ((dt/(density * 2 * dx)) * (p[i+1][j] - p[i-1][j])) +
            (viscosity * ( ((dt/(dx*dx)) * (un[i+1][j] - ((2*un[i][j]) + un[i-1][j]))) +
                           ((dt/(dy*dy)) * (un[i][j+1] - ((2*un[i][j]) + un[i][j-1])))))
        );

        v[i][j] = (
          (vn[i][j]) -
          (un[i][j] * (dt/dx) * (vn[i][j] - vn[i-1][j])) -
          (vn[i][j] * (dt/dy) * (vn[i][j] - vn[i][j-1])) -
          ((dt/(density * 2 * dx)) * (p[i][j+1] - p[i][j-1])) +
          (viscosity * ( ((dt/(dx*dx)) * (vn[i+1][j] - ((2*vn[i][j]) + vn[i-1][j]))) +
                         ((dt/(dy*dy)) * (vn[i][j+1] - ((2*vn[i][j]) + vn[i][j-1])))))
        );
      }
    }

    /*  Boundary Conditions  */
    for(int i = 0; i < NX; i++)
    {
      u[i][0] = 0;
      u[i][NY-1] = 1;     // should be lid vel
      v[i][0] = 0;
      v[i][NY-1] = 0;
    }
    for(int i = 0; i < NY; i++)
    {
      u[0][i] = 0;
      u[NX-1][i] = 0;
      v[0][i] = 0;
      v[NX-1][i] = 0;
    }
  }

}


void pressure_rhs(double **b, double **u, double **v, double rho, double dt,
  double dx, double dy)
{
  // ignore boundaries
  double inv_dt = 1/dt;
  for(int i = 1; i < (NX-1); i++)
  {
    for(int j = 1; j < (NY-1); j++)
    {
      // could square lines 2 and 4
      // could calculate some of these values and store them every loop (would save
      //    some calculations [2-3 per loop])
      b[i][j] = rho * (
        (inv_dt * ( ((u[i+1][j] - u[i-1][j]) / (2*dx)) + ((v[i][j+1] - v[i][j-1]) / (2*dy)) ) ) -
        ( ((u[i+1][j] - u[i-1][j]) / (2*dx)) * ((u[i+1][j] - u[i-1][j]) / (2*dx)) ) -
        (2 * ((u[i][j+1] - u[i][j-1])/(2*dy)) * ((v[i+1][j] - v[i-1][j])/(2*dx)) ) -
        ( ((v[i][j+1] - v[i][j-1])/(2*dy)) * ((v[i][j+1] - v[i][j-1])/(2*dy)) )
      );
    }
  }
}

void pressure_poisson(double **p, double dx, double dy, double** b, int max_iter, double eps)
{
  // epsilon convergence
  // do {
  //
  // } while()
  double **pn = alloc_table(NX, NY);
  double dy2 = dy*dy;
  double dx2 = dx*dx;
  for(int it = 0; it < max_iter; it++)
  {
    dynamic_memcpy(p, pn, NX, NY);
    for(int i = 1; i < (NX-1); i++)
    {
      for(int j = 1; j < (NY-1); j++)
      {
        p[i][j] = (
          ((((pn[i+1][j] + pn[i-1][j]) * dy2) + ((pn[i][j+1] + pn[i][j-1]) * dx2)) / (2*(dx2 + dy2))) -
          ((dx2*dy2)/(2*(dx2 + dy2))) * b[i][j]
        );
      }
    }
  }

  for(int i = 0; i < NX; i++)
  {
    p[NX-1][i] = p[NX-2][i];    // dp/dx = 0 @ x = x_max
    p[0][i] = p[1][i];          // dp/dx = 0 @ x = 0
  }
  for(int i = 0; i < NY; i++)
  {
    p[i][0] = p[i][1];          // dp/dy = 0 @ y =0
    p[i][NY-1] = 0;             // p = 0 @ y = 2
  }
}

// copies a dynamically allcoated 2D array (a) to another dynamically allocated
// 2D array
void dynamic_memcpy(double **a, double **b, int rows, int cols)
{
  for(int i = 0; i < rows; i++)
  {
    memcpy(b[i], a[i], cols*sizeof(double));
  }
}

double** alloc_table(int numx, int numy)
{
  double **tab = new double*[numx];
  for(int i = 0; i < numx; i++)
  {
      tab[i] = new double[numy];
  }

  for(int i = 0; i < numx; i++)
  {
    for(int j = 0; j < numy; j++)
    {
      tab[i][j] = 0;
    }
  }

  return tab;
}
