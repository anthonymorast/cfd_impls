class parameters
{
public:
    void compile();

    /* simulation variables */
    int nx = 0;                 // grid spaces x-direction
    int ny = 0;                 // grid spaces y-direction
    int nt = 0;                 // number of time steps
    double dx;
    double dy;
    double dt = 0;

    /* Pressure Convergence */
    int max_iter_pp = 0;        // max iterations to try for convergence of Pressure Poisson
    double pp_epsilon = 0;     // stops Pressure Poisson convergence

    /* geometry variables */
    double x_max = 0;           // rhs for x
    double y_max = 0;           // rhs for y
    double cylinder_x = -1;     // cylinder center x
    double cylinder_y = -1;     // cylinder center y
    double cylinder_radius = -1;// radius of cylinder

    /* Re = (length*density*velocity)/viscosity */
    double density = 0;
    double viscosity = 0;

    /* Boundary Condition function pointers; */
    /* each takes velocity and pressure fields as input as well as the grid dims */
    void (*top)(double **, double**, double **, double **, double **, int, int);    // top b.c.
    void (*rhs)(double **, double**, double **, double **, double **, int, int);    // right-hand side b.c.
    void (*bot)(double **, double**, double **, double **, double **, int, int);    // bottom b.c.
    void (*lhs)(double **, double**, double **, double **, double **, int, int);    // left-hand side b.c.
};

/*
*
*   The compile function is called after the simulation parameters are specified.
*   The function will set other variables used by the simulation. The variables
*   will be public so setting them before calling compile will do nothing but
*   setting them afterwards will be catastrophic (probably).
*
*/
void parameters::compile()
{
    dx = x_max / (nx-1);
    dy = y_max / (ny-1);
}
