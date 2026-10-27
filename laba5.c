#include <stdio.h>
#include <math.h>

double f(double t) {
    if (fabs(t) < 1e-10) {
        return 0;
    }
    return (cos(t) - 1.0) / t;
}


double method(double a, double b, int n) {   
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 3 == 0) {
            sum += 2.0 * f(x);
        } 
        else {
            sum += 3.0 * f(x);
        }
    }
    
    return (3.0 * h * sum) / 8.0;
}

int main() {
    double eps = 1e-9;
    double up;
    if (scanf("%lf", &up) != 1) return 1;

    double maskerony = 0.5772156649015328;
    double lnx = log(up);

    int n = 3;
    double result1, result2;

    do {
        result1 = maskerony + lnx + method(0, up, n);
        
        n += 3;
        
        result2 = maskerony + lnx + method(0, up, n);
        
    } while (fabs(result2 - result1) > eps);

    printf("%.15lf\n", result2);

    return 0;
}