
#include <stdio.h>
#include <math.h>
#include "mat.h"
#include <stdlib.h>

#define R 6378
#define N 3602
#define GM 398600.5
#define alt 230

double rad(double degree)
{
     return (degree * M_PI) / 180.0; // prevod na radiány
}

double Rtotal(double altitude) // spocitanie vysky
{
     return (R + alt);
}

// vypocet koordinatov v jednom cykle
void calculate_coordinates(MAT *coordinatesS, MAT *coordinatesX, MAT *coordinatesE, double *B, double *L, int n, double altitude)
{
     for (int i = 0; i < n; i++)
     {
          double Brad = rad(B[i]);
          double Lrad = rad(L[i]);

          // S koordinaty
          ELEM(coordinatesS, i, 0) = Rtotal(0) * cos(Brad) * cos(Lrad); // ELEM je element matice na prislusnom mieste
          ELEM(coordinatesS, i, 1) = Rtotal(0) * cos(Brad) * sin(Lrad);
          ELEM(coordinatesS, i, 2) = Rtotal(0) * sin(Brad);

          // X koordinaty
          ELEM(coordinatesX, i, 0) = Rtotal(altitude) * cos(Brad) * cos(Lrad);
          ELEM(coordinatesX, i, 1) = Rtotal(altitude) * cos(Brad) * sin(Lrad);
          ELEM(coordinatesX, i, 2) = Rtotal(altitude) * sin(Brad);

          // E koordinaty
          ELEM(coordinatesE, i, 0) = cos(Brad) * cos(Lrad);
          ELEM(coordinatesE, i, 1) = cos(Brad) * sin(Lrad);
          ELEM(coordinatesE, i, 2) = sin(Brad);
     }
}

double distance(double *x1, double *x2)
{
     return sqrt(pow(x2[0] - x1[0], 2) + pow(x2[1] - x1[1], 2) + pow(x2[2] - x1[2], 2));
}

void calculate_distance_matrix(MAT *distanceMatrix, MAT *coordinatesX, MAT *coordinatesS, int n)
{
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               double p1[3] = {ELEM(coordinatesX, i, 0), ELEM(coordinatesX, i, 1), ELEM(coordinatesX, i, 2)};
               double p2[3] = {ELEM(coordinatesS, j, 0), ELEM(coordinatesS, j, 1), ELEM(coordinatesS, j, 2)};
               ELEM(distanceMatrix, i, j) = distance(p1, p2);
          }
     }
}

void calculate_Aij(MAT *A, MAT *distanceMatrix, MAT *distanceVectors, MAT *coordinatesE, int n)
{
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               double rij = ELEM(distanceMatrix, i, j);
               double dotProduct = 0.0;
               dotProduct = rij * ELEM(coordinatesE, i, j);
               ELEM(A, i, j) = (1 / pow(rij, 3)) - ((3 * pow(dotProduct, 2)) / pow(rij, 5));
          }
     }
}

void load_data(const char *filename, MAT *B, MAT *L, MAT *H, MAT *dg, MAT *f, int n)
{
     FILE *file = fopen(filename, "r");
     if (file == NULL)
     {
          printf("Error opening file.\n");
          exit(1);
     }

     for (int i = 0; i < n; i++)
     {
          double bi, li, hi, dgi, fi;
          fscanf(file, "%lf %lf %lf %lf %lf", &bi, &li, &hi, &dgi, &fi);
          ELEM(B, i, 0) = bi;
          ELEM(L, i, 0) = li;
          ELEM(H, i, 0) = hi;
          ELEM(dg, i, 0) = dgi;
          ELEM(f, i, 0) = fi;
     }

     fclose(file);
}

// Matrix multiplication function
// pridavam novu funkciu
void mat_multiply(MAT *result, const MAT *A, const MAT *B)
{
    if (A->cols != B->rows)
    {
        printf("Matrix dimensions do not match for multiplication.\n");
        return;
    }

    // Initialize the result matrix to zero
    mat_zero(result);

    for (unsigned int i = 0; i < A->rows; i++)
    {
        for (unsigned int j = 0; j < B->cols; j++)
        {
            for (unsigned int k = 0; k < A->cols; k++)
            {
                ELEM(result, i, j) += ELEM(A, i, k) * ELEM(B, k, j);
            }
        }
    }
}


int main()
{
     int n = 3602;
     MAT *coordinatesS = mat_create_with_type(n, 3);
     MAT *coordinatesX = mat_create_with_type(n, 3);
     MAT *coordinatesE = mat_create_with_type(n, 3);
     MAT *distanceMatrix = mat_create_with_type(n, n);
     MAT *distanceVectors = mat_create_with_type(n, n * 3);
     MAT *A = mat_create_with_type(n, n);

     MAT *B = mat_create_with_type(n, 1);  // Matice pre šírku (B)
     MAT *L = mat_create_with_type(n, 1);  // Matice pre dĺžku (L)
     MAT *H = mat_create_with_type(n, 1);  // Matice pre výšku (H)
     MAT *dg = mat_create_with_type(n, 1); // Matice pre dg
     MAT *f = mat_create_with_type(n, 1);  // Matice pre f

     load_data("C:\\Users\\puvak\\OneDrive\\Počítač\\Timovy projekt\\GOCE-4\\BL-3602.dat", B, L, H, dg, f, n); // Úprava názvu súboru podľa potreby

     double B_vals[n], L_vals[n];
     for (int i = 0; i < n; i++) 
     {
        B_vals[i] = ELEM(B, i, 0);  // Get B values
        L_vals[i] = ELEM(L, i, 0);  // Get L values
     }

     // Výpočet súradníc
     calculate_coordinates(coordinatesS, coordinatesX, coordinatesE,B_vals, L_vals, n, alt);

     // Výpočet matice vzdialeností
     calculate_distance_matrix(distanceMatrix, coordinatesX, coordinatesS, n);

     // Výpočet matice Aij
     calculate_Aij(A, distanceMatrix, distanceVectors, coordinatesE, n);

     // Riešenie pre alpha
     MAT *alpha = mat_create_with_type(n, 1);
     MAT *dgM = mat_create_with_type(n, 1);

     // Fill dgM with data (you can use load_data or fill it manually)
    for (int i = 0; i < n; i++) {
        ELEM(dgM, i, 0) = ELEM(dg, i, 0); // Assuming dgM should hold dg values
    }

    // Solve the system A * alpha = dgM
    if (mat_division(A, dgM, alpha) == SUCCESS) {
        printf("Alpha vector:\n");
        mat_print(alpha); // Print the result
    } else {
        printf("Failed to solve the system.\n");
    }

    // Calculate the final vector u (this assumes some operation with alpha)
    MAT *u = mat_create_with_type(n, 1);
    // Assuming u = A * alpha for demonstration (replace with your actual calculation)
    mat_multiply(u, A, alpha); // If you want to compute u = A * alpha

    //Print the result
    printf("Final vector u:\n");
    mat_print(u);

     MAT *u1 = mat_create_with_type(n, 1);
    mat_division(u1,A,alpha);
    printf("Final vector u:\n");
    mat_print(u);



    // Free allocated memory
    mat_destroy(coordinatesS);
    mat_destroy(coordinatesX);
    mat_destroy(coordinatesE);
    mat_destroy(distanceMatrix);
    mat_destroy(distanceVectors);
    mat_destroy(A);
    mat_destroy(alpha);
    mat_destroy(dgM);
    mat_destroy(u);
    mat_destroy(B);
    mat_destroy(L);
    mat_destroy(H);
    mat_destroy(dg);
    mat_destroy(f);

    return 0;
}