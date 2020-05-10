#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "the2.h"


/*
INPUT:
    int *row_count: vertical size of the matrix to be assigned (passed by reference from the main function)
    int *row_size: horizontal size of the matrix to be assigned (passed by reference from the main function)

OUTPUT:
    double **matrix: created data matrix

METHOD:
    This function creates the matrix streamed via stdin. 
    During this process, assigns row_count and row_size variables.
*/
double **initialize_the_data(int *row_count, int *row_size) {
    int i, j;
    double **matrix;
    double data;
    matrix = (double**) malloc(1 * sizeof(double*));
    for (i = 0;; i++){
        matrix = (double**) realloc(matrix, (i+2) * sizeof(double*));
        matrix[i] = (double*) malloc(1 * sizeof(double));
        for (j = 0;; j++){
            scanf("%lf",&data);
            if (data == -1){
                break;
            }
            else {
                matrix[i] = (double*) realloc(matrix[i], (j+2) * sizeof(double));
                matrix[i][j] = data;
                *row_size = j+1;
            }
        }
        if (j == 0){
            *row_count = i;
            break;
        }
    }
    return matrix; 
    
}

/*
INPUT:
    double **matrix: data matrix
    int n: number of rows to be printed
    int row_size: horizontal size of the matrix

METHOD:
    This function prints first n row of the matrix.
*/
void print_first_n_row(double **matrix, int n, int row_size) {
    int i, j;
    for (i = 0; i < n; i++){
        for (j = 0; j <= row_size; j++){
            if ( j == row_size && i == n-1){
                break;
            }
            else if (j == row_size){
                printf("\n");
            }
            else if (j == row_size-1){
                printf("%.4lf",matrix[i][j]);
            }
            else{
                printf("%.4lf ",matrix[i][j]);
            }
            
        }
    }
}

/*
INPUT:
    double **matrix: data matrix
    int row_size: horizontal size of the data matrix
    int row1: index of the first row in the calculation
    int row2: index of the second row in the calculation

METHOD:
    This function calculates dot product of the row1 and the row2 and prints it in the following format:
        Dot product of rows <row1>, <row2>: <dot_product>
*/
void calculate_dot_product(double **matrix, int row_size, int row1, int row2) {
    int i;
    double result = 0;
    for (i = 0;i < row_size; i++){
        result = result + (matrix[row1][i] * matrix[row2][i]);
    }
    printf("Dot product of rows %d, %d: %.4lf",row1,row2,result);
}

/*
INPUT:
    double **matrix: data matrix
    int row_count: vertical size of the data matrix
    int row_size: horizontal size of the data matrix

OUTPUT:
    double **covariance_matrix: Created covariance matrix with size of row_size X row_size

METHOD:
    This function creates covariance matrix of the original data matrix and returns it.

*/
float mean(int row, double **matrix, int row_count){
    int i;
    float result = 0;
    for (i = 0; i < row_count; i++){
        result += matrix[i][row];
    }
    return result / row_count;
}
float calculate(int row1, int row2, double **matrix, int row_count, int row_size){
    int i;
    float mean1, mean2;
    float result = 0;
    for (i = 0; i < row_count; i++){
        result += (matrix[i][row1] - mean(row1, matrix, row_count)) * (matrix[i][row2] - mean(row2, matrix, row_count));
    }
    return result / (row_count - 1);
}

double **calculate_covariance_matrix(double **matrix, int row_count, int row_size) { 
    int i, j;
    double **cov_matrix;
    cov_matrix = (double**) malloc(row_size * sizeof(double*));
    for (i = 0; i < row_size; i++){
        cov_matrix[i] = (double*) malloc(row_size * sizeof(double));
        for (j = 0; j < row_size; j++){
            cov_matrix[i][j] = calculate(i, j, matrix, row_count, row_size);
        }
    }
    return cov_matrix;

    
}

/*
INPUT:
    double **matrix: data matrix
    int row_count: vertical size of the data matrix
    int row_size: horizontal size of the data matrix

OUTPUT:
    double **result: Created result matrix with size of row_size X row_size

METHOD:
    This function calculates x^T * x. First finds the transpose of the original data matrix and apply matrix multiplication.
    At the end it returns the calculated matrix with size of row_size X row_size.

*/

double **calculate_x_transpose_times_x(double **matrix, int row_count, int row_size) { 
    double **result_matrix;
    int i, j, k;
    double result;
    result_matrix = (double**) malloc(row_size * sizeof(double*));
    for (i = 0; i < row_size; i++){
        result_matrix[i] = (double*) malloc(row_size * sizeof(double));
        for (j = 0; j < row_size; j++){
            result = 0;
            for (k = 0; k < row_count; k++){
                result += matrix[k][i] * matrix[k][j];
            }
            result_matrix[i][j] = result;
            
            
        }
    }
    return result_matrix;
    
}

/*
INPUT:
    double **matrix: data matrix
    int *group_count: number of the groups to be assigned (passed by reference from the main function)
    int row_count: vertical size of the data matrix
    int row_size: horizontal size of the data matrix
    int group_column: index of the column to apply grouping over
    int operation: index of the operation to apply during grouping
        SUM -> 0
        AVG -> 1
        MAX -> 2
        MIN -> 3

OUTPUT:
    double **grouped_matrix: Created grouped matrix with size of group_count X row_size

METHOD:
    This function applies grouping over desired column index, creates a grouped matrix. It returns this grouped matrix.
    During this process it calculates group count and assigns to the variable.

*/
int in(double data, double *list, int array_size){
    int i;
    for (i = 0; i < array_size; i++){
        if (list[i] == data){
            return 1;
        }
    }
    return 0;
}
double **group_by(double **matrix, int *group_count, int row_count, int row_size, int group_column, int operation) { 
    double **grouped_matrix;
    int i, j, k, l, t, n, m;
    int counter = 0;
    double data;
    double **deletematrix;
    double *used;
    int minic = 0;
    used = (double*) malloc(row_count * sizeof(double));
    grouped_matrix = (double**) malloc(sizeof(double*));
    for (i = 0; i < row_count; i++){
        data = matrix[i][group_column];
        if (i != 0 && in(data, used, row_count) == 1 ){
            continue;
        }
        used[counter] = data;
        grouped_matrix = (double**) realloc(grouped_matrix, (counter+1) * sizeof(double*));
        grouped_matrix[counter] = (double*) malloc(row_size * sizeof(double));

        for (k = 0; k < row_size; k++){
            if (k == group_column){
                grouped_matrix[counter][k] = data;
            }
            else{
                grouped_matrix[counter][k] = 0;
            }
        }
        for (j = 0; j < row_count; j++){
            if (matrix[j][group_column] == data){
                minic++;
                if (operation == 0){
                    for (l = 0; l < row_size; l++){
                        if ( l != group_column){
                            grouped_matrix[counter][l] += matrix[j][l];
                        }
                    }
                    
                    
                }
                if (operation == 1){
                    for (l = 0; l < row_size; l++){
                        if ( l != group_column){
                            grouped_matrix[counter][l] += matrix[j][l];
                        }
                    }
                        
                }
                if (operation == 2){
                    for (n = 0; n < row_size; n++){
                        if (n != group_column){
                            if (grouped_matrix[counter][n] == 0){
                                grouped_matrix[counter][n] = matrix[j][n];
                            }
                            if (grouped_matrix[counter][n] < matrix[j][n]){
                                grouped_matrix[counter][n] = matrix[j][n];
                            }
                        }
                    }
                    
                }
                if (operation == 3){
                    for (m = 0; m < row_size; m++){
                        if (n != group_column){
                            if (grouped_matrix[counter][m] == 0){
                                grouped_matrix[counter][m] = matrix[j][m];
                            }
                            if (grouped_matrix[counter][m] > matrix[j][m]){
                                grouped_matrix[counter][m] = matrix[j][m];
                            }
                        }
                    }
                    
                }
                
            }
        }
        if (operation == 1){
            for (t = 0; t < row_size; t++){
                if (t != group_column){
                    grouped_matrix[counter][t] /= minic;
                }
            }
        }
        minic = 0;
        counter++;
    }
    *group_count = counter;
    return grouped_matrix;
    
}

/*
INPUT:
    double **matrix: data matrix
    int row_count: vertical size of the data matrix
    int row_size: horizontal size of the data matrix
    double **kernel: kernel matrix
    int kernel_height: vertical size of the kernel matrix
    int kernel_width: horizontal size of the kernel matrix

OUTPUT:
    double **convoluted_matrix: Created convoluted matrix

METHOD:
    This function applies convolution over data matrix using kernel matrix and returns created matrix.

*/
double **convolution(double **matrix, int row_count, int row_size, double **kernel, int kernel_height, int kernel_width) { 
    double **convoluted;
    int new_height = row_count - kernel_height + 1;
    int new_width = row_size - kernel_width + 1;
    int i, j ,k ,l ,n;
    double result;
    convoluted = (double**) malloc(new_height * sizeof(double*));
    for (i = 0; i < new_height; i++){
        convoluted[i] = (double*) malloc(new_width * sizeof(double));
    }
    for (j = 0; j < new_height; j++){
        for (k = 0; k < new_width; k++){
            for (l = 0; l < kernel_height; l++){
                for (n = 0; n < kernel_width; n++){
                    result += kernel[l][n] * matrix[j+l][k+n];
                }
            }
            convoluted[j][k] = result;
            result = 0;
            
        }
    }
    return convoluted;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}
