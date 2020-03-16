##
# Example 4.1 in section 4.3 from H K Versteeg and W Malalasekera's book.
##
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    ## Given
    nx = 100
    k = 1000
    T_A = 100
    T_B = 500
    A = 10*10e-3
    xmin = 0
    xmax = 0.5

    ## Derived from given
    dx = (xmax - xmin)/(nx)     # note distance from T_A to x_1 = dx/2 similar for T_B to x_n
                                # this makes dx 0.1 rather than 0.125. This is also
                                # demonstrated in the illustration wherein the faces
                                # are half the distance between nodal values.
    x = np.linspace(xmin+dx, xmax, nx)      # actual x coords
    T = np.zeros(nx)

    aw = (k/dx)*A
    ae = aw
    ap = aw + ae
    udiff = 1
    while udiff > 1e-50:
        Tn = T.copy()
        T[0] = ((k*A/dx)*Tn[1] + (2*k*A/dx)*T_A)/((k*A/dx) + ((2*k*A)/dx))
        T[-1] = ((k*A/dx)*Tn[-2] + ((2*k*A/dx)*T_B))/((k*A/dx) + (2*k*A/dx))
        T[1:-1] = (aw*Tn[0:-2] + ae*Tn[2:])/ap
        udiff = (np.sum(T) - np.sum(Tn)) / np.sum(T)
    print(T)
    plt.plot(x, T)
    plt.show()
