#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(){
    FILE *file = fopen("insert_value.txt", "r");
    
    if (file == NULL){
        printf("Error opening file!\n");
    }

    int size;
    fscanf(file, "%d", &size);

    double x[size];
    double y[size];
    for (int i = 0; i < size; i++){
        fscanf(file, "%lf", &x[i]);
    }

    for (int i = 0; i < size; i++){
        fscanf(file, "%lf", &y[i]);
    }

    for (int i = 0; i < size; i++){
        printf("%.2lf ", x[i]);
        printf("%.2lf ", y[i]);
        printf("\n");
    }
    
    printf("\n");
    
    double x_1 = 0, x_2 = 0, x_3 = 0, x_4 = 0, x_5 = 0, x_6 = 0, x_2y = 0, x_1y = 0, x_3y = 0, y_1 = 0;

    for (int i = 0; i < size; i++){
        x_2 += pow(x[i], 2);
        x_3 += pow(x[i], 3);
        x_4 += pow(x[i], 4);
        x_5 += pow(x[i], 5);
        x_6 += pow(x[i], 6);
        x_1 += x[i];
    }

    for (int i = 0; i < size; i++){
        x_2y += pow(x[i], 2) * y[i];
        x_1y += x[i] * y[i];
        x_3y += pow(x[i], 3) * y[i];
        y_1 += y[i];
    }
    
    int n = 4;
    int m = 5;

    double matrix[n][m]; 
    matrix[0][0] = size;
    matrix[0][1] = x_1;
    matrix[0][2] = x_2;
    matrix[0][3] = x_3;
    matrix[0][4] = y_1;
    matrix[1][0] = x_1;
    matrix[1][1] = x_2;
    matrix[1][2] = x_3;
    matrix[1][3] = x_4;
    matrix[1][4] = x_1y;
    matrix[2][0] = x_2;
    matrix[2][1] = x_3;
    matrix[2][2] = x_4;
    matrix[2][3] = x_5;
    matrix[2][4] = x_2y;
    matrix[3][0] = x_3;
    matrix[3][1] = x_4;
    matrix[3][2] = x_5;
    matrix[3][3] = x_6;
    matrix[3][4] = x_3y;

    double matrix2[n][m];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            matrix2[i][j] = matrix[i][j];
        }
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            printf("%.2lf ", matrix2[i][j]);
        }
        printf("\n");
    }


    printf("\n");
 
    for (int k = 0; k < n; k++){
        int max_row = k;
        double max_val = fabs(matrix2[k][k]);
        for (int i = k + 1; i < n; i++){
            if (fabs(matrix2[i][k]) > max_val){
                max_val = fabs(matrix2[i][k]);
                max_row = i;
            }
        }

        if (max_row != k){
            for (int j = 0; j < m; j++){
                double temp = matrix2[k][j];
                matrix2[k][j] = matrix2[max_row][j];
                matrix2[max_row][j] = temp;
            }
        }

        if (fabs(matrix2[k][k]) < 1e-12){
            fprintf(stderr, "Error: Matrix is singular or nearly singular!\n");
            return 1;
        }

        // Исключение элементов ниже главной диагонали
        for (int i = k + 1; i < n; i++){
            double factor = matrix2[i][k] / matrix2[k][k];
            for (int j = k; j < m; j++){
                matrix2[i][j] -= factor * matrix2[k][j];
            }
        }
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            printf("%.2lf ", matrix2[i][j]);
        }
        printf("\n");
    }

    double *result = (double *)malloc(n * sizeof(double));
    if (result == NULL){
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }

    for (int i = n - 1; i >= 0; i--){
        result[i] = matrix2[i][n];
        for (int j = i + 1; j < n; j++){
            result[i] -= matrix2[i][j] * result[j];
        }
        result[i] /= matrix2[i][i];
    }

    printf("\n");

    for (int i = 0; i < n; i++){
        printf("a[%d] = %.4lf\n", i, result[i]);
    }

    double a0, a1, a2, a3;
    a0 = result[0];
    a1 = result[1];
    a2 = result[2];
    a3 = result[3];

    double y_new[size];
    FILE *file2 = fopen("output.txt", "w");
    if (file2 == NULL){
        fprintf(stderr, "Error: Failed to open output.txt\n");
        return 1;
    }

    for (int i = 0; i < size; i++){
        y_new[i] = a0 + a1 * x[i] + a2 * x[i] * x[i] + a3 * x[i] * x[i] * x[i];
    }

    double error[size];
    for (int i = 0; i < size; i++){
        error[i] = (y_new[i] - y[i]);
    }

    fprintf(file2, "Start function sqrt(x)\n");
    fprintf(file2, "polinom: 0.5060+0.5331x-0.0452x^2+0.0018x^3\n");
    fprintf(file2, "x\t\t y\t\t y_new\t error\n");
    for (int i = 0; i < size; i++){
        fprintf(file2, "%.4lf\t %.4lf\t %.4lf\t %.4lf\n", x[i], y[i], y_new[i], error[i]);
    }


    fclose(file);

    return 0;

}

