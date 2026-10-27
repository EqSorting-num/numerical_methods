#include <stdio.h>
#include <math.h>

double f(double x) {
    return sqrt(1 + x * x) + exp(-2 * x); 
}

int main() {
    double a = 0.0;         
    double b = 1.0;         
    double eps = 0.001;   

    int n = 10;

    while ((b-a) > eps){
        double h = (b-a)/n;
        double min_x = a;
        double min_f = f(a);
        int min_index = 0;

        for (int i = 0; i <= n; i++){
            double current_x = a + i*h;
            double current_f = f(current_x);

            if (current_f < min_f){
                min_x = current_x;
                min_f = current_f;
                min_index = i;
            }   
        }

        double new_a = min_x - h;
        double new_b = min_x + h;

        if (new_a < a){
            new_a = a;
        }

        if (new_b > b){
            new_b = b;
        }

        a = new_a;
        b = new_b;

        if (h <= 0){
            break;
        }
    }

    double final_x = (a+b)/2;
    double final_f = f(final_x);
    printf("min_x = %.4f\n", final_x);
    printf("min_f = %.4f\n", final_f);
    return 0;


}
