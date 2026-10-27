
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.001

static const double B1 = 100.0, B2 = 0.0, T1 = 0.0, T2 = 0.0;

static double K(double x) {
    return x * x + 1.0; 
}
static double q(double x) { 
    return x + 2.0; 
}
static double f(double x) { 
    return cos(x); 
}

static void solve(int N, double *y)
{
    double h = 1.0 / N;
    double *A = malloc((N + 1) * sizeof(double));
    double *C = malloc((N + 1) * sizeof(double));
    double *Bc = malloc((N + 1) * sizeof(double));
    double *D = malloc((N + 1) * sizeof(double));
    double *al = malloc((N + 1) * sizeof(double));
    double *be = malloc((N + 1) * sizeof(double));
    int i;


    {
        double a1 = K(0.5 * h);
        C[0]  = a1 / h + B1 + 0.5 * h * q(0.0);
        Bc[0] = -a1 / h;
        A[0]  = 0.0;
        D[0]  = T1 + 0.5 * h * f(0.0);
    }

    for (i = 1; i < N; i++) {
        double x  = i * h;
        double ai = K(x - 0.5 * h);
        double a1 = K(x + 0.5 * h);
        A[i]  = -ai / (h * h);
        C[i]  = (ai + a1) / (h * h) + q(x);
        Bc[i] = -a1 / (h * h);
        D[i]  = f(x);
    }

    {
        double aN = K(1.0 - 0.5 * h);
        A[N]  = -aN / h;
        C[N]  = aN / h + B2 + 0.5 * h * q(1.0);
        Bc[N] = 0.0;
        D[N]  = T2 + 0.5 * h * f(1.0);
    }

    al[0] = -Bc[0] / C[0];
    be[0] = D[0] / C[0];
    for (i = 1; i <= N; i++) {
        double den = C[i] + A[i] * al[i - 1];
        al[i] = -Bc[i] / den;
        be[i] = (D[i] - A[i] * be[i - 1]) / den;
    }
    y[N] = be[N];
    for (i = N - 1; i >= 0; i--)
        y[i] = al[i] * y[i + 1] + be[i];

    free(A); free(C); free(Bc); free(D); free(al); free(be);
}

int main(void)
{
    int N = 10;
    const int NMAX = 1 << 22;
    double err = 1e9, prev_err = 0.0;

    printf("Grid method, eps = %g\n", EPS);
    printf("%8s %12s %10s\n", "N", "grade err", "attitude");

    while (N <= NMAX) {
        int N2 = 2 * N, i;
        double *y1 = malloc((N + 1) * sizeof(double));
        double *y2 = malloc((N2 + 1) * sizeof(double));

        solve(N, y1);
        solve(N2, y2);

        err = 0.0;
        for (i = 0; i <= N; i++) {
            double d = fabs(y1[i] - y2[2 * i]) / 3.0;
            if (d > err) err = d;
        }
        if (prev_err > 0.0)
            printf("%8d %12.3e %10.2f\n", N, err, prev_err / err);
        else
            printf("%8d %12.3e %10s\n", N, err, "-");
        prev_err = err;

        if (err < EPS) {
            double h = 1.0 / N2;
            int step = N2 / 10;
            if (step < 1) step = 1;

            printf("\nstep h = %g (N = %d)\n\n", h, N2);
            printf("%8s %14s\n", "x", "U(x)");
            for (i = 0; i <= N2; i += step)
                printf("%8.3f %14.8f\n", i * h, y2[i]);
            printf("\nU(0) = %.8f, U(1) = %.8f\n", y2[0], y2[N2]);
            free(y1); free(y2);
            return 0;
        }
        free(y1); free(y2);
        N = N2;
    }
    printf("The specified accuracy has not been achieved\n");
    return 0;
}


