#include <stdio.h>
#include <math.h>

#define N 25

double load[N];

/* Linear interpolation */
double f(double t)
{
    int i = (int)t;

    if(i < 0) i = 0;
    if(i >= N-1) return load[N-1];

    double frac = t - i;

    return load[i] + frac * (load[i+1] - load[i]);
}

/* Gauss-Legendre Quadrature */
double gaussQuadrature(int nPoints)
{
    int bins = 24;           // fixed as required
    double a = 0.0;
    double b = 24.0;
    double h = (b-a)/bins;

    double total = 0.0;

    double x[5], w[5];

    switch(nPoints)
    {
        case 1:
            x[0] = 0.0;
            w[0] = 2.0;
            break;

        case 2:
            x[0] = -0.5773502692;
            x[1] =  0.5773502692;

            w[0] = 1.0;
            w[1] = 1.0;
            break;

        case 3:
            x[0] = -0.7745966692;
            x[1] =  0.0;
            x[2] =  0.7745966692;

            w[0] = 0.5555555556;
            w[1] = 0.8888888889;
            w[2] = 0.5555555556;
            break;

        case 4:
            x[0] = -0.8611363116;
            x[1] = -0.3399810436;
            x[2] =  0.3399810436;
            x[3] =  0.8611363116;

            w[0] = 0.3478548451;
            w[1] = 0.6521451549;
            w[2] = 0.6521451549;
            w[3] = 0.3478548451;
            break;

        case 5:
            x[0] = -0.9061798459;
            x[1] = -0.5384693101;
            x[2] =  0.0;
            x[3] =  0.5384693101;
            x[4] =  0.9061798459;

            w[0] = 0.2369268851;
            w[1] = 0.4786286705;
            w[2] = 0.5688888889;
            w[3] = 0.4786286705;
            w[4] = 0.2369268851;
            break;

        default:
            return -1;
    }

    for(int i=0;i<bins;i++)
    {
        double xl = i;
        double xr = i+1;

        double mid  = (xl+xr)/2.0;
        double half = (xr-xl)/2.0;

        double sum = 0.0;

        for(int j=0;j<nPoints;j++)
        {
            double t = mid + half*x[j];
            sum += w[j]*f(t);
        }

        total += half*sum;
    }

    return total;
}

int main()
{
    printf("Input 25 load values:\n");

    for(int i=0;i<25;i++)
        scanf("%lf",&load[i]);

    printf("\n=== Gauss n-Point Comparison ===\n\n");

    for(int n=1;n<=5;n++)
    {
        printf("Gauss %d-Point : %.2lf MWh\n",
               n,
               gaussQuadrature(n));
    }

    return 0;
}