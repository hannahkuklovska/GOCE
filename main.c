
#include <stdio.h>
#include <math.h>
#include "mat.h"

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
               for (int k = 0; k < 3; k++)
               {
                    dotProduct += ELEM(distanceVectors, i, j, k) * ELEM(coordinatesE, i, k);
               }
               ELEM(A, i, j) = (1 / pow(rij, 3)) - ((3 * pow(dotProduct, 2)) / pow(rij, 5));
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
     MAT *distanceVectors = mat_create_with_type(n, n, 3);
     MAT *A = mat_create_with_type(n, n);

     // Arrays for B, L (read from file)
     double B[n], L[n];
     // Load data into B, L arrays

     // Coordinate calculations
     calculate_coordinates(coordinatesS, coordinatesX, coordinatesE, B, L, n, alt);

     // Distance matrix calculation
     calculate_distance_matrix(distanceMatrix, coordinatesX, coordinatesS, n);

     // Aij calculation
     calculate_Aij(A, distanceMatrix, distanceVectors, coordinatesE, n);

     // Solve for alpha
     MAT *alpha = mat_create_with_type(n, 1);
     MAT *dgM = mat_create_with_type(n, 1);
     // Fill dgM array with values
     // Solve A * alpha = dgM

     // Final u vector calculation
     MAT *u = mat_create_with_type(n, 1);
     // Calculate u

     // Print the result
     mat_print(u);

     // Free memory
     mat_destroy(coordinatesS);
     mat_destroy(coordinatesX);
     mat_destroy(coordinatesE);
     mat_destroy(distanceMatrix);
     mat_destroy(distanceVectors);
     mat_destroy(A);
     mat_destroy(alpha);
     mat_destroy(dgM);
     mat_destroy(u);

     return 0;
}