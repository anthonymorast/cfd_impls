/*
    The arrays here can be a little confusing; at least they were when I implemented
    this solver.

    I designed the primitive variable arrays to hold the data as though it were an
    actual grid overlaying a cartesian coordinate system. That is, the first row
    of the arrays are the values at the "top" of the boundary. The last row, i.e.
    row NY-1, is the "bottom" boundary. The same is true for the X directions:
    the 0th index is the LHS of the boundary and the NX-1 index is the RHS.

    Conceptualizing an array stored in memory we would get something like this:

    _ _ _ _ _ _ _ _ ... _     row 1
    _ _ _ _ _ _ _ _ ... _     row 2
    _ _ _ _ _ _ _ _ ... _     row 3
    etc.

    So the first index U[0][:] is the first Row or the "first" Y coordinate which,
    as described above is the Y-MAX = NY row. The same is true for the second index,
    U[:][0] which would be the LHS of the boundary. Thus, U[0][0] is the top-left
    point in the grid, U[NY-1][NX-1] is the bottom right point in the grid, etc.
*/
#include "bitmap_image.hpp"
#include "utils.h"
#include "parameters.h"

#include <math.h>
#include <time.h>

struct flow_data
{
    double **u;
    double **v;
    double **un;
    double **vn;
    double **p;
    double **pn;
    double **b;
    double **F;
};

void top(double** p, double** u, double** v, double **b, double **F, int nx, int ny);
void bot(double** p, double** u, double** v, double **b, double **F, int nx, int ny);
void rhs(double** p, double** u, double** v, double **b, double **F, int nx, int ny);
void lhs(double** p, double** u, double** v, double **b, double **F, int nx, int ny);

void apply_boundaries(parameters &p, flow_data &d);
parameters create_params();
parameters get_params_from_file(string filename);
void iterate(parameters &p, flow_data &d);
int converge_pressure(parameters &p, flow_data &d);
void build_poisson_rhs(parameters &p, flow_data &d);
void solve_momentum(parameters &p, flow_data &d);
void create_image(bitmap_image img, double **data, string filename);
bool check_nan(double **data, int ny, int nx);

int main(int argc, char *argv[])
{
    parameters p = create_params();
    flow_data data;

    data.u = alloc_table(p.ny, p.nx);
    data.v = alloc_table(p.ny, p.nx);
    data.un = alloc_table(p.ny, p.nx);
    data.vn = alloc_table(p.ny, p.nx);
    data.p = alloc_table(p.ny, p.nx);
    data.pn = alloc_table(p.ny, p.nx);
    data.b = alloc_table(p.ny, p.nx);
    data.F = alloc_table(p.ny, p.nx);

    iterate(p, data);
}

void iterate(parameters &p, flow_data &d)
{
    bitmap_image image(p.nx, p.ny);
    int total_iterations_pp = 0;
    int image_save_iter = int(0.25/p.dt);   // every .25 seconds = every image_save_iter iterations
    clock_t start, end;
	  double timeTaken = 0;

    for(int t = 0; t < p.nt; t++)
    {
        cout << "Iteration " << t << " of " << p.nt << "..." << endl;
        dynamic_memcpy(d.u, d.un, p.ny, p.nx);
        dynamic_memcpy(d.v, d.vn, p.ny, p.nx);

        start = clock();
        build_poisson_rhs(p, d);
        end = clock();
        timeTaken = double(end-start) / double(CLOCKS_PER_SEC);
        cout << "\tBuild PPE RHS: " << timeTaken << " sec" << endl;

        start = clock();
        total_iterations_pp += converge_pressure(p, d);
        end = clock();
        timeTaken = double(end-start) / double(CLOCKS_PER_SEC);
        cout << "\tConverge/Iterate Pressure: " << timeTaken << " sec" << endl;

        start = clock();
        solve_momentum(p, d);
        end = clock();
        timeTaken = double(end-start) / double(CLOCKS_PER_SEC);
         cout << "\tSolve Momentum Equations: " << timeTaken << " sec" << endl;

        start = clock();
        apply_boundaries(p, d);
        end = clock();
        timeTaken = double(end-start) / double(CLOCKS_PER_SEC);
        cout << "\tApply Boundaries: " << timeTaken << " sec" << endl;

        if(check_nan(d.u, p.ny, p.nx) ||
           check_nan(d.v, p.ny, p.nx) ||
           check_nan(d.p, p.ny, p.nx))
            return;

        if(t % image_save_iter == 0)
        {
            cout << t << " iterations complete of "
                 << p.nt << ", writing results to file..." << endl;
            string filename = "iter_" + to_string(t);
            create_image(image, d.u, "imgs/u/" + filename + "_u.bmp");
            create_image(image, d.v, "imgs/v/" + filename + "_v.bmp");
            create_image(image, d.p, "imgs/p/" + filename + "_p.bmp");

            array_to_csv(d.u, p.nx, p.ny, "data/u/" + filename + "_u.csv", "u", t);
            array_to_csv(d.v, p.nx, p.ny, "data/v/" + filename + "_v.csv", "v", t);
            array_to_csv(d.p, p.nx, p.ny, "data/p/" + filename + "_p.csv", "p", t);
        }
    }

    cout << "Average PP Iterations: " << ((double)total_iterations_pp)/((double)p.nt) << endl;
}

void create_image(bitmap_image img, double** data, string filename)
{
    double max = -999999999;
    double min = 999999999;

    for(int i = 0; i < img.height(); i++)
    {
        for(int j = 0; j < img.width(); j++)
        {
            if(data[i][j] < min)
            {
                min = data[i][j];
            }
            if(data[i][j] > max)
            {
                max = data[i][j];
            }
        }
    }

    // http://www.andrewnoske.com/wiki/Code_-_heatmaps_and_color_gradients
    int NUM_COLORS = 4;
    float color[NUM_COLORS][3] = { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
    int idx1, idx2;
    float fracBetween = 0;


    img.clear();
    for(int j = 0; j < img.height(); j++)
    {
        for(int i = 0; i < img.width(); i++)
        {
            double value = 0;
            if((max - min) != 0)
            {
                value = (data[j][i] - min / (max - min));
            }

            if(value <= 0)
            {
                idx1 = idx2 = 0;
            }
            else if (value >= 1)
            {
                idx1 = idx2 = NUM_COLORS-1;
            }
            else
            {
                value = value * (NUM_COLORS - 1);
                idx1 = floor(value);
                idx2 = idx1 + 1;
                fracBetween = value - float(idx1);
            }
            float red = (color[idx2][0] - color[idx1][0])*fracBetween + color[idx1][0];
            float green = (color[idx2][1] - color[idx1][1])*fracBetween + color[idx1][1];
            float blue = (color[idx2][2] - color[idx1][2])*fracBetween + color[idx1][2];
            red *= 255; green *= 255; blue *= 255;
            img.set_pixel(i, j, red, green, blue);
        }
    }
    img.save_image(filename);
}

bool check_nan(double** data, int ny, int nx)
{
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            if(isnan(data[j][i]))
            {
                cout << "NAN found in array!" << endl;
                cout << "\tindex: (" << j << "," << i << ")" << endl;
                return true;
            }
        }
    }
    return false;
}

void solve_momentum(parameters &p, flow_data &d)
{
    double au = p.dt/(p.density * 2 * p.dx);
    double av = p.dt/(p.density * 2 * p.dy);
    double c = p.dt/(p.dx*p.dx);
    double e = p.dt/(p.dy*p.dy);

    double dtdx = p.dt/p.dx;
    double dtdy = p.dt/p.dy;

    for(int j = 1; j < (p.ny-1); j++)
    {
        for(int i = 1; i < (p.nx-1); i++)
        {
            d.u[j][i] = (
                d.un[j][i] -
                (d.un[j][i] * dtdx * (d.un[j][i] - d.un[j][i-1])) -
                (d.vn[j][i] * dtdy * (d.un[j][i] - d.un[j-1][i])) -
                (au * (d.p[j][i+1] - d.p[j][i-1])) +
                (p.viscosity * ((c * (d.un[j][i+1] - 2*d.un[j][i] + d.un[j][i-1]))
                              + (e * (d.un[j+1][i] - 2*d.un[j][i] + d.un[j-1][i])))) +
                (d.F[j][i]*p.dt)
            );

            d.v[j][i] = (
                d.vn[j][i] -
                (d.un[j][i] * dtdx * (d.vn[j][i] - d.vn[j][i-1])) -
                (d.vn[j][i] * dtdy * (d.vn[j][i] - d.vn[j-1][i])) -
                (av * (d.p[j+1][i] - d.p[j-1][i])) +
                (p.viscosity * ((c * (d.vn[j][i+1] - 2*d.vn[j][i] + d.vn[j][i-1]))
                              + (e * (d.vn[j+1][i] - 2*d.vn[j][i] + d.vn[j-1][i]))))
            );
        }
    }
}

/**
* Returns the number of iterations it took to converge the PP
**/
int converge_pressure(parameters &p, flow_data &d)
{
    int iter = 0;
    double max_change;      // max change on the grid
    double dx2 = p.dx*p.dx;
    double dy2 = p.dy*p.dy;

    // setting the pressure field to 1 before convergence took avg iters from
    // 218.685 to 4635.92

    // double average = 0;
    // for(int i = 0; i < p.ny; i++)
    // {
    //     for(int j = 0; j < p.nx; j++)
    //     {
    //         average += d.p[i][j];
    //     }
    // }
    // average /= (p.ny * p.nx);
    //
    // // add the average pressure to the pressure field
    // // adding (average * .0001) to the pressure field takes the average Iterations
    // // up to 386.639
    // for(int i = 0; i < p.ny; i++)
    // {
    //     for(int j = 0; j < p.nx; j++)
    //     {
    //         // d.p[i][j] += (average*.000001);
    //     }
    // }
    do
    {
        max_change = 0;
        dynamic_memcpy(d.p, d.pn, p.ny, p.nx);
        for(int j = 1; j < (p.ny-1); j++)
        {
            for(int i = 1; i < (p.nx-1); i++)
            {
                d.p[j][i] = (
                    (((dy2*(d.pn[j][i+1] + d.pn[j][i-1]))+(dx2*(d.pn[j+1][i] + d.pn[j-1][i])))/(2*(dx2 + dy2))) -
                    ((dx2*dy2)/(2*(dx2 + dy2))) *
                    d.b[j][i]
                );
                max_change = max(fabs(d.p[j][i] - d.pn[j][i]), max_change);
            }
        }
        iter++;

        // need to apply pressure boundaries, others shouldn't matter
        apply_boundaries(p, d);
    } while(iter < p.max_iter_pp && max_change > p.pp_epsilon);

    return iter;
}

void build_poisson_rhs(parameters &p, flow_data &d)
{
    // save that precious compute time
    double idt = 1/p.dt;
    double dxdx = 2*p.dx;
    double dydy = 2*p.dy;

    for(int j = 1; j < (p.ny-1); j++)
    {
        for(int i = 1; i < (p.nx-1); i++)
        {
            d.b[j][i] = p.density * (
                (idt*(((d.u[j][i+1] - d.u[j][i-1])/(dxdx)) + ((d.v[j+1][i] - d.v[j-1][i])/(dydy)))) -
                pow((d.u[j][i+1] - d.u[j][i-1])/(dxdx), 2) -
                (2*(((d.u[j+1][i] - d.u[j-1][i])/(dydy)) * ((d.v[j][i+1] - d.v[j][i-1])/(dxdx)))) -
                pow((d.v[j+1][i] - d.v[j-1][i])/(dydy),2)
            );
        }
    }
}

parameters create_params()
{
    parameters p;

    p.nx = 501;
    p.ny = 501;
    p.nt = 50000;
    p.dt = .00001;
    p.max_iter_pp = 500;
    p.pp_epsilon = 0; //.00000001;
    p.x_max = 2;
    p.y_max = 2;
    p.density = 1.0;
    p.viscosity = 0.1;
    p.top = top;
    p.rhs = rhs;
    p.bot = bot;
    p.lhs = lhs;

    p.compile();

    return p;
}


parameters get_params_from_file(string filename)
{
    parameters p;
    return p;
}

void apply_boundaries(parameters &p, flow_data &d)
{
    p.rhs(d.p, d.u, d.v, d.b, d.F, p.nx, p.ny);
    p.lhs(d.p, d.u, d.v, d.b, d.F, p.nx, p.ny);
    p.bot(d.p, d.u, d.v, d.b, d.F, p.nx, p.ny);
    p.top(d.p, d.u, d.v, d.b, d.F, p.nx, p.ny);
}

void top(double** p, double** u, double** v, double **b, double **F, int nx, int ny)
{
    for(int i = 0; i < nx; i++)
    {
        p[0][i] = 0;                // p = 0 @ y = y_max
        u[0][i] = 1;                // u = 1 at lid
        v[0][i] = 0;
    }
}

void bot(double** p, double** u, double** v, double **b, double **F, int nx, int ny)
{
    for(int i = 0; i < nx; i++)
    {
        p[ny-1][i] = p[ny-2][i];    // dp/dy = 0 @ y = 0
        u[ny-1][i] = 0;
        v[ny-1][i] = 0;
    }
}

void rhs(double** p, double** u, double** v, double **b, double **F, int nx, int ny)
{
    for(int i = 0; i < ny; i++)
    {
        p[i][nx-1] = p[i][nx-2];
        u[i][nx-1] = 0;
        v[i][nx-1] = 0;
    }
}

void lhs(double** p, double** u, double** v, double **b, double **F, int nx, int ny)
{
    for(int i = 0; i < ny; i++)
    {
        p[i][0] = p[i][1];
        u[i][0] = 0;
        v[i][0] = 0;
    }
}
