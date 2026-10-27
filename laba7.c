#include <stdio.h>
#include <math.h>

#define EPS 0.001
#define ALPHA 0.1  // фиксированный размер шага

double f(double x1, double x2) {
    return 2.0 * (x1 - 1.0) * (x1 - 1.0) + x2 * x2 + 2.0 * x1 * x2;
}

void compute_gradient(double x1, double x2, double *g1, double *g2) {
    *g1 = 4.0 * x1 + 2.0 * x2 - 4.0;
    *g2 = 2.0 * x1 + 2.0 * x2;
}

int main() {
    double x1 = 0.0;
    double x2 = 0.0;
    
    double g1, g2;
    int iter = 0;

    while (1) {
        compute_gradient(x1, x2, &g1, &g2);

        double grad_norm = sqrt(g1 * g1 + g2 * g2);

        if (grad_norm < EPS) {
            printf("Сходимость достигнута!\n");
            break;
        }

        // Делаем шаг с фиксированным размером ALPHA
        x1 = x1 - ALPHA * g1;
        x2 = x2 - ALPHA * g2;

        iter++;
        
        // Ограничение на число итераций
        if (iter > 10000) {
            printf("Превышено максимальное число итераций.\n");
            break;
        }
    }

    printf("Решение: x1 = %.5f, x2 = %.5f, f = %.5f\n", x1, x2, f(x1, x2));
    printf("Итераций: %d\n", iter);

    return 0;
}