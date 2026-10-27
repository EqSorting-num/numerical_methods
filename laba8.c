
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A    1.0
#define B    2.0
#define U0   0.8
#define DU0  0.6
#define EPS  0.000001

static void f(double x, const double y[2], double out[2])
{
    out[0] = y[1];
    out[1] = 2.0 * y[0] / (x * x);
}

static void rk4_step(double x, const double y[2], double h, double yn[2])
{
    double k1[2], k2[2], k3[2], k4[2], t[2];
    int i;

    f(x, y, k1);
    for (i = 0; i < 2; i++) t[i] = y[i] + 0.5 * h * k1[i];
    f(x + 0.5 * h, t, k2);
    for (i = 0; i < 2; i++) t[i] = y[i] + 0.5 * h * k2[i];
    f(x + 0.5 * h, t, k3);
    for (i = 0; i < 2; i++) t[i] = y[i] + h * k3[i];
    f(x + h, t, k4);
    for (i = 0; i < 2; i++)
        yn[i] = y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}


static void hamming(int N, double (*y)[2], double (*fv)[2],double (*p)[2], double (*c)[2])
{
    double h = (B - A) / N;
    int n, i;

    y[0][0] = U0;
    y[0][1] = DU0;

    for (n = 0; n < 3; n++)
        rk4_step(A + n * h, y[n], h, y[n + 1]);

    for (n = 0; n <= 3; n++)
        f(A + n * h, y[n], fv[n]);

    for (i = 0; i < 2; i++) {
        p[3][i] = y[3][i];
        c[3][i] = y[3][i];
    }

    for (n = 3; n < N; n++) {
        double xn1 = A + (n + 1) * h;
        double m[2], fm[2];

        for (i = 0; i < 2; i++) {
            p[n + 1][i] = y[n - 3][i]
                        + 4.0 * h / 3.0 * (2.0 * fv[n][i] - fv[n - 1][i] + 2.0 * fv[n - 2][i]);
            m[i] = p[n + 1][i] + 112.0 / 121.0 * (c[n][i] - p[n][i]);
        }

        f(xn1, m, fm);

        for (i = 0; i < 2; i++) {
            c[n + 1][i] = (9.0 * y[n][i] - y[n - 2][i]
                           + 3.0 * h * (fm[i] + 2.0 * fv[n][i] - fv[n - 1][i])) / 8.0;
            y[n + 1][i] = c[n + 1][i] + 9.0 / 121.0 * (p[n + 1][i] - c[n + 1][i]);
        }

        f(xn1, y[n + 1], fv[n + 1]);
    }
}

static double exact(double x)
{
    return 7.0 / 15.0 * x * x + 1.0 / (3.0 * x);
}

int main(void)
{
    int N = 4;               
    const int NMAX = 1 << 20;
    double (*y1)[2], (*f1)[2], (*p1)[2], (*c1)[2];
    double (*y2)[2], (*f2)[2], (*p2)[2], (*c2)[2];
    double err = 1e9;
    int i;

    y1 = f1 = p1 = c1 = NULL;
    y2 = f2 = p2 = c2 = NULL;

    printf("Hamming method, eps = %g\n", EPS);
    printf("%8s %12s\n", "N", "grade err");

    while (N <= NMAX) {
        int N2 = 2 * N;

        y1 = realloc(y1, (N + 1) * sizeof *y1);
        f1 = realloc(f1, (N + 1) * sizeof *f1);
        p1 = realloc(p1, (N + 1) * sizeof *p1);
        c1 = realloc(c1, (N + 1) * sizeof *c1);
        y2 = realloc(y2, (N2 + 1) * sizeof *y2);
        f2 = realloc(f2, (N2 + 1) * sizeof *f2);
        p2 = realloc(p2, (N2 + 1) * sizeof *p2);
        c2 = realloc(c2, (N2 + 1) * sizeof *c2);

        hamming(N,  y1, f1, p1, c1);
        hamming(N2, y2, f2, p2, c2);

        err = 0.0;
        for (i = 0; i <= N; i++) {
            double d = fabs(y1[i][0] - y2[2 * i][0]) / 15.0;
            if (d > err) err = d;
        }
        printf("%8d %12.3e\n", N, err);

        if (err < EPS) {
            double h = (B - A) / N2;
            double maxerr = 0.0;

            printf("\nЧNumber of steps: %d, h = %g\n\n", N2, h);
            printf("%10s %14s %14s %14s\n", "x", "u(x)", "u'(x)", "|u - u_точн|");

            int step = N2 / 10;
            if (step < 1) step = 1;
            for (i = 0; i <= N2; i++) {
                double x = A + i * h;
                double e = fabs(y2[i][0] - exact(x));
                if (e > maxerr) maxerr = e;
                if (i % step == 0 || i == N2)
                    printf("%10.5f %14.8f %14.8f %14.3e\n", x, y2[i][0], y2[i][1], e);
            }
            printf("\nu(2) = %.8f (exact %.8f)\n", y2[N2][0], exact(B));
            printf("Actully error: %.3e\n", maxerr);
            break;
        }
        N = N2;
    }

    if (err >= EPS)
        printf("The specified accuracy has not been achieved\n");

    free(y1); free(f1); free(p1); free(c1);
    free(y2); free(f2); free(p2); free(c2);
    return 0;
}