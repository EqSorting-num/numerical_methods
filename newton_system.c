#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TOLERANCE 1e-4
#define MAX_ITER 100

/* Structure for 2x2 matrix */
typedef struct {
    double a[2][2];
} Matrix2x2;

/* Structure for 2D vector */
typedef struct {
    double v[2];
} Vector2;

/* System of equations F(x,y) = 0 */
Vector2 f(double x, double y) {
    Vector2 F;
    F.v[0] = sin(x + y) - 1.2 * x - 0.2;
    F.v[1] = x*x + y*y - 1;
    return F;
}

/* Jacobian matrix J */
Matrix2x2 jacobian(double x, double y) {
    Matrix2x2 J;
    J.a[0][0] = cos(x + y) - 1.2;
    J.a[0][1] = cos(x + y);
    J.a[1][0] = 2 * x;
    J.a[1][1] = 2 * y;
    return J;
}

/* Determinant of 2x2 matrix */
double det(Matrix2x2 M) {
    return M.a[0][0] * M.a[1][1] - M.a[0][1] * M.a[1][0];
}

/* Inverse of 2x2 matrix */
Matrix2x2 inverse(Matrix2x2 M) {
    Matrix2x2 inv;
    double d = det(M);
    double inv_d = 1.0 / d;
    
    inv.a[0][0] =  M.a[1][1] * inv_d;
    inv.a[0][1] = -M.a[0][1] * inv_d;
    inv.a[1][0] = -M.a[1][0] * inv_d;
    inv.a[1][1] =  M.a[0][0] * inv_d;
    
    return inv;
}

/* Multiply matrix by vector */
Vector2 multiply(Matrix2x2 M, Vector2 v) {
    Vector2 result;
    result.v[0] = M.a[0][0] * v.v[0] + M.a[0][1] * v.v[1];
    result.v[1] = M.a[1][0] * v.v[0] + M.a[1][1] * v.v[1];
    return result;
}

/* Vector norm (Euclidean) */
double norm(Vector2 v) {
    return sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]);
}

/* Newton's method for systems */
int newton_method(double x0, double y0, double *x, double *y) {
    double x_curr = x0, y_curr = y0;
    double x_next, y_next;
    int iter = 0;
    
    printf("\nSystem:\n");
    printf("  sin(x + y) - 1.2x = 0.2\n");
    printf("  x^2 + y^2 = 1\n");
    printf("\nTarget accuracy: %.0e\n", TOLERANCE);
    printf("Initial guess: x0 = %.1f, y0 = %.1f\n", x0, y0);
    printf("------------------------------------------------------------\n");
    printf("%5s %12s %12s %12s %12s\n", "Iter", "x", "y", "||F||", "||d||");
    printf("------------------------------------------------------------\n");
    
    for (iter = 0; iter < MAX_ITER; iter++) {
        Vector2 F = f(x_curr, y_curr);
        
        Matrix2x2 J = jacobian(x_curr, y_curr);
        
        double d = det(J);
        if (fabs(d) < 1e-15) {
            printf(" Jacobian is singular at iteration %d\n", iter);
            return -1;
        }
        
        Matrix2x2 J_inv = inverse(J);
        Vector2 negF;
        negF.v[0] = -F.v[0];
        negF.v[1] = -F.v[1];
        Vector2 delta = multiply(J_inv, negF);
        
        x_next = x_curr + delta.v[0];
        y_next = y_curr + delta.v[1];
        
        double norm_F = norm(F);
        double norm_delta
        ;
        
        printf("%5d %12.8f %12.8f %12.8f %12.8f\n", 
               iter + 1, x_curr, y_curr, norm_F, norm_delta);
        
        if (norm_F < TOLERANCE && norm_delta < TOLERANCE) {
            printf("\nSolution: x = %.8f, y = %.8f\n", x_next, y_next);
            *x = x_next;
            *y = y_next;
            return iter + 1;
        }
        
        x_curr = x_next;
        y_curr = y_next;
    }
    

    printf("\nMaximum iterations reached without full convergence\n");
    printf("  Last approximation: x = %.8f, y = %.8f\n", x_curr, y_curr);
    *x = x_curr;
    *y = y_curr;
    return iter;
}

/* Verify the solution */
void verify_solution(double x, double y) {
    
    double eq1 = sin(x + y) - 1.2 * x;
    double eq2 = x*x + y*y;
    
    printf("  Equation 1: sin(x + y) - 1.2x = %.8f (expected: 0.2)\n", eq1);
    printf("  Equation 2: x^2 + y^2 = %.8f (expected: 1)\n", eq2);
    
    double error1 = fabs(eq1 - 0.2);
    double error2 = fabs(eq2 - 1);
    
    printf("\n  Error in eq1: |%.8f - 0.2| = %.8e\n", eq1, error1);
    printf("  Error in eq2: |%.8f - 1| = %.8e\n", eq2, error2);
    
}

int main() {
    double x1, y1, x2, y2;
    int iters1, iters2;
    
    printf("\nInitial guess (0.5, 0.5)\n");
    iters1 = newton_method(0.5, 0.5, &x1, &y1);
    verify_solution(x1, y1);
    
    printf("\nInitial guess (-0.5, 0.5)\n");
    iters2 = newton_method(-0.5, 0.5, &x2, &y2);
    verify_solution(x2, y2);
    
    printf("%20s %20s %15s\n", "Initial Guess", "Solution (x, y)", "Iterations");
    printf("%20s (%8.8f, %8.8f) %15d\n", "(0.5, 0.5)", x1, y1, iters1);
    printf("%20s (%8.8f, %8.8f) %15d\n", "(-0.5, 0.5)", x2, y2, iters2);
    
    return 0;
}
