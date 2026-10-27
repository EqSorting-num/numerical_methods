#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void mat_vec_mul(double **A, double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = 0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
}

double vec_norm(double *x, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}

void vec_normalize(double *x, int n) {
    double norm = vec_norm(x, n);
    for (int i = 0; i < n; i++) {
        x[i] /= norm;
    }
}

int solve_system(double **A, double *b, double *y, int n) {
    double **a = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        a[i] = (double *)malloc((n + 1) * sizeof(double));
        for (int j = 0; j < n; j++) {
            a[i][j] = A[i][j];
        }
        a[i][n] = b[i];
    }

    for (int k = 0; k < n; k++) {
        int max_row = k;
        double max_val = fabs(a[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(a[i][k]) > max_val) {
                max_val = fabs(a[i][k]);
                max_row = i;
            }
        }

        if (max_val < 1e-12) {
            printf("The matrix is ​​singular.\n");
            for (int i = 0; i < n; i++) free(a[i]);
            free(a);
            return 0;
        }

        if (max_row != k) {
            double *temp = a[k];
            a[k] = a[max_row];
            a[max_row] = temp;
        }

        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k] / a[k][k];
            for (int j = k; j <= n; j++) {
                a[i][j] -= factor * a[k][j];
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        y[i] = a[i][n];
        for (int j = i + 1; j < n; j++) {
            y[i] -= a[i][j] * y[j];
        }
        y[i] /= a[i][i];
    }

    for (int i = 0; i < n; i++) free(a[i]);
    free(a);
    return 1;
}

void power_method(double **A, int n, double eps, double *lambda_max, double *vector) {
    double *z = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        z[i] = 1.0;
    }
    vec_normalize(z, n);

    double lambda_old = 0;
    double lambda_new = 0;
    
    do {
        lambda_old = lambda_new; 
        
        mat_vec_mul(A, z, y, n);
        lambda_new = vec_norm(y, n);
        
        vec_normalize(y, n);
        for (int i = 0; i < n; i++) {
            z[i] = y[i];
        }
        
    } while (fabs(lambda_new - lambda_old) > eps);

    *lambda_max = lambda_new;
    for (int i = 0; i < n; i++) {
        vector[i] = z[i];
    }

    free(z);
    free(y);
}

void inverse_iteration(double **A, int n, double eps, double *lambda_min, double *vector) {
    double *z = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));
    double *z_old = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        z[i] = 1.0;
    }
    vec_normalize(z, n);

    double lambda_old = 0;
    double lambda_new = 0;
    int iter = 0;

    do {
        lambda_old = lambda_new;

        for (int i = 0; i < n; i++) {
            z_old[i] = z[i];
        }

        if (!solve_system(A, z, y, n)) {
            printf("Ошибка решения системы!\n");
            return;
        }

        double norm_z = vec_norm(z_old, n);
        double norm_y = vec_norm(y, n);
        lambda_new = norm_z / norm_y;

        vec_normalize(y, n);
        for (int i = 0; i < n; i++) {
            z[i] = y[i];
        }

    } while (fabs(lambda_new - lambda_old) > eps && iter < 1000);

    *lambda_min = lambda_new;
    for (int i = 0; i < n; i++) {
        vector[i] = z[i];
    }

    free(z);
    free(y);
    free(z_old);
}

int main() {
    double eps = 1e-3;

    FILE *file = fopen("matrix.txt", "r");
    if (file == NULL) {
        printf("Error file open\n");
        return 1;
    }

    int n;
    fscanf(file, "%d", &n);

    double **A = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%lf", &A[i][j]);
        }
    }
    fclose(file);

    double lambda_max, lambda_min;

    double *vec_max = (double *)malloc(n * sizeof(double));
    double *vec_min = (double *)malloc(n * sizeof(double));

    power_method(A, n, eps, &lambda_max, vec_max);
    printf("\nmax lambda = %.6lf\n", lambda_max);
    printf("eigenvector: [");
    for (int i = 0; i < n; i++) {
        printf("%.6lf", vec_max[i]);
        if (i < n - 1) printf(", ");
    }
    printf("]\n");

    inverse_iteration(A, n, eps, &lambda_min, vec_min);
    printf("\nmin lambda = %.6lf\n", lambda_min);
    printf("eigenvector: [");
    for (int i = 0; i < n; i++) {
        printf("%.6lf", vec_min[i]);
        if (i < n - 1) printf(", ");
    }
    printf("]\n");

    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(vec_max);
    free(vec_min);

    return 0;
}