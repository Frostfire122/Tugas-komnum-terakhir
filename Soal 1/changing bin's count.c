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

/* Left Riemann */
double leftRiemann(int bins)
{
    double a = 0.0;
    double b = 24.0;
    double h = (b-a)/bins;

    double sum = 0;

    for(int i=0;i<bins;i++)
        sum += f(a + i*h);

    return sum*h;
}

/* Right Riemann */
double rightRiemann(int bins)
{
    double a = 0.0;
    double b = 24.0;
    double h = (b-a)/bins;

    double sum = 0;

    for(int i=1;i<=bins;i++)
        sum += f(a + i*h);

    return sum*h;
}

/* Midpoint */
double midpointRiemann(int bins)
{
    double a = 0.0;
    double b = 24.0;
    double h = (b-a)/bins;

    double sum = 0;

    for(int i=0;i<bins;i++)
        sum += f(a + (i+0.5)*h);

    return sum*h;
}

/* Trapezoidal */
double trapezoidal(int bins)
{
    double a = 0.0;
    double b = 24.0;
    double h = (b-a)/bins;

    double sum = (f(a)+f(b))/2.0;

    for(int i=1;i<bins;i++)
        sum += f(a+i*h);

    return sum*h;
}

/* Gauss 2-Point */
double gauss2(int bins)
{
    double a = 0.0;
    double b = 24.0;
    double h = (b-a)/bins;

    double total = 0;

    for(int i=0;i<bins;i++)
    {
        double xl = a + i*h;
        double xr = xl + h;

        double mid  = (xl+xr)/2.0;
        double half = (xr-xl)/2.0;

        double x1 = mid - half/sqrt(3.0);
        double x2 = mid + half/sqrt(3.0);

        total += half*(f(x1)+f(x2));
    }

    return total;
}

int main()
{
    printf("Input 25 hourly load values:\n");

    for(int i=0;i<25;i++)
        scanf("%lf",&load[i]);

    int bins[] = {1,2,4,6,8,12,24};
    int totalBins = 7;

    printf("\n=============================================\n");
    printf("Energy Estimation (MWh)\n");
    printf("=============================================\n");

    for(int i=0;i<totalBins;i++)
    {
        int b = bins[i];

        printf("\nBins = %d\n", b);
        printf("Left Riemann  : %.2lf\n", leftRiemann(b));
        printf("Right Riemann : %.2lf\n", rightRiemann(b));
        printf("Midpoint      : %.2lf\n", midpointRiemann(b));
        printf("Trapezoidal   : %.2lf\n", trapezoidal(b));
        printf("Gauss 2-Point : %.2lf\n", gauss2(b));
    }

    return 0;
}