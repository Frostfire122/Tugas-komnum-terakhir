#include <stdio.h>
#include <math.h>

#define N 25

/* =====================================================
   RIEMANN SUM
   ===================================================== */
double riemannSum(double y[], int n, double h)
{
    double sum = 0.0;

    for(int i = 0; i < n - 1; i++)
        sum += y[i];

    return sum * h;
}

/* =====================================================
   TRAPEZOIDAL RULE
   ===================================================== */
double trapezoidalRule(double y[], int n, double h)
{
    double sum = y[0] + y[n - 1];

    for(int i = 1; i < n - 1; i++)
        sum += 2.0 * y[i];

    return (h / 2.0) * sum;
}

/* =====================================================
   SIMPSON'S RULE
   ===================================================== */
double simpsonRule(double y[], int n, double h)
{
    int intervals = n - 1;
    double sum = y[0] + y[n - 1];

    for(int i = 1; i < intervals; i++)
    {
        if(i % 2 == 0)
            sum += 2.0 * y[i];
        else
            sum += 4.0 * y[i];
    }

    return (h / 3.0) * sum;
}

/* =====================================================
   LINEAR INTERPOLATION
   ===================================================== */
double interpolate(double x,
                   double x0,
                   double x1,
                   double y0,
                   double y1)
{
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

/* =====================================================
   GAUSSIAN QUADRATURE (2-POINT)
   ===================================================== */
double gaussianQuadrature(double y[], int intervals)
{
    double total = 0.0;

    for(int i = 0; i < intervals; i++)
    {
        double a = i;
        double b = i + 1;

        double x1 = (a + b) / 2.0 - (b - a) / (2.0 * sqrt(3.0));
        double x2 = (a + b) / 2.0 + (b - a) / (2.0 * sqrt(3.0));

        double f1 = interpolate(x1, a, b, y[i], y[i + 1]);
        double f2 = interpolate(x2, a, b, y[i], y[i + 1]);

        total += ((b - a) / 2.0) * (f1 + f2);
    }

    return total;
}

/* =====================================================
   COMBINED METHOD (RTSG)
   ===================================================== */
double RTSG_Method(double r,
                   double t,
                   double s,
                   double g)
{
    return 0.10 * r
         + 0.20 * t
         + 0.40 * s
         + 0.30 * g;
}

/* =====================================================
   MAIN
   ===================================================== */
int main()
{
    double load[N];

    for(int i = 0; i < N; i++)
    {
        scanf("%lf", &load[i]);
    }

    double h = 1.0;
    int intervals = N - 1;

    double riemann   = riemannSum(load, N, h);
    double trapezoid = trapezoidalRule(load, N, h);
    double simpson   = simpsonRule(load, N, h);
    double gaussian  = gaussianQuadrature(load, intervals);

    double combined = RTSG_Method(
        riemann,
        trapezoid,
        simpson,
        gaussian
    );

    printf("Riemann Sum : %.2f\n", riemann);
    printf("Trapezoidal Rule : %.2f\n", trapezoid);
    printf("Simpson's Rule : %.2f\n", simpson);
    printf("Gaussian Quadrature : %.2f\n", gaussian);
    printf("RTSG Method (Riemann + Trapezoidal + Simpson + Gaussian) : %.2f\n", combined);

    return 0;
}