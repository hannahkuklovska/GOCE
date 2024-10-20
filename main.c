
#include <stdio.h>
#include <math.h>
#include "mat.h"
#include <stdlib.h>

#define R 6378
#define N 3602
#define GM 398600.5
#define alt 230
#define dGM 0.0000270141

MAT *mat_create_with_type(unsigned int rows, unsigned int cols)
{
     MAT *mat = (MAT *)malloc(sizeof(MAT));

     if (mat == NULL)
     {
          return NULL;
     }
     mat->rows = rows;
     mat->cols = cols;
     mat->elem = (float *)malloc(rows * cols * sizeof(float));
     if (mat->elem == NULL)
     {
          free(mat);
          return NULL;
     }
     return mat;
}

MAT *mat_zero(MAT *null)
{
     unsigned i, j;

     for (i = 0; i < null->rows; i++)
     {
          for (j = 0; j < null->cols; j++)
          {
               ELEM(null, i, j) = 0.0;
          }
     }
     return null;
}

MAT *mat_create_by_file(char *filename)
{
     unsigned int rows, cols;
     char prve[2];

     FILE *file = fopen(filename, "rb");
     if (file == NULL)
     {
          return NULL;
     }

     fread(prve, sizeof(char), 2, file);
     if (prve[0] != 'M' || prve[1] != '1')
     {
          fclose(file);
          return NULL;
     }

     fread(&rows, sizeof(unsigned int), 1, file);
     fread(&cols, sizeof(unsigned int), 1, file);

     // vytvorenie matice
     MAT *mat = mat_create_with_type(rows, cols);
     if (mat == NULL)
     {
          fclose(file);
          return NULL;
     }

     // prečítanie matice
     fread(mat->elem, sizeof(float), rows * cols, file);
     fclose(file);
     return mat;
}

char mat_save(MAT *mat, char *filename)
{
     FILE *file = fopen(filename, "wb");
     if (file == NULL)
     {
          return FAILURE;
     }
     fwrite("M1", sizeof(char), 2, file);
     fwrite(&(mat->rows), sizeof(unsigned int), 1, file);
     fwrite(&(mat->cols), sizeof(unsigned int), 1, file);
     fwrite(mat->elem, sizeof(float), mat->rows * mat->cols, file);
     fclose(file);
     return SUCCESS;
}

void mat_destroy(MAT *mat)
{
     if (mat != NULL)
     {
          free(mat->elem);
          free(mat);
     }
}

void mat_unit(MAT *mat)
{
     unsigned int i, j;

     for (i = 0; i < mat->rows; i++)
     {
          for (j = 0; j < mat->cols; j++)
          {
               if (i == j)
               {
                    ELEM(mat, i, j) = 1.0;
               }
               else
               {
                    ELEM(mat, i, j) = 0.0;
               }
          }
     }
}

void mat_random(MAT *mat)
{
     unsigned int i, j;

     for (i = 0; i < mat->rows; i++)
     {
          for (j = 0; j < mat->cols; j++)
          {
               ELEM(mat, i, j) = -1.0 + 2 * ((float)rand()) / RAND_MAX;
          }
     }
}

void mat_print(MAT *mat)
{
     unsigned int i, j;

     for (i = 0; i < mat->rows; i++)
     {
          for (j = 0; j < mat->cols; j++)
          {
               printf("%.30f\t", ELEM(mat, i, j));
          }
          printf("\n");
     }
}

MAT *mat_invert(MAT *input_matrix)
{
     unsigned int dimension;
     MAT *augmented_matrix;
     MAT *inverse_matrix;
     int i, j, k;
     float pivot_value, factor;

     if (input_matrix->rows != input_matrix->cols)
     {
          return NULL;
     }

     dimension = input_matrix->rows;
     augmented_matrix = mat_create_with_type(dimension, 2 * dimension);
     inverse_matrix = mat_create_with_type(dimension, dimension);

     if (inverse_matrix == NULL || augmented_matrix == NULL)
     {
          mat_destroy(augmented_matrix);
          mat_destroy(inverse_matrix);
          return NULL;
     }

     for (i = 0; i < dimension; i++)
     {
          for (j = 0; j < dimension; j++)
          {
               ELEM(augmented_matrix, i, j) = ELEM(input_matrix, i, j);

               if (j == i)
               {
                    ELEM(augmented_matrix, i, j + dimension) = 1.0;
               }
               else
               {
                    ELEM(augmented_matrix, i, j + dimension) = 0.0;
               }
          }
     }
     // Gaussova eliminacna metoda
     for (i = 0; i < dimension; i++)
     {
          // najdenie pivota
          pivot_value = ELEM(augmented_matrix, i, i);
          if (pivot_value == 0)
          {
               mat_destroy(augmented_matrix);
               mat_destroy(inverse_matrix);
               return NULL;
          }
          // normalizacia pivota
          for (j = 0; j < 2 * dimension; j++)
          {
               ELEM(augmented_matrix, i, j) /= pivot_value;
          }

          for (k = 0; k < dimension; k++)
          {
               if (k == i)
               {
                    continue;
               }
               factor = ELEM(augmented_matrix, k, i);
               for (j = 0; j < 2 * dimension; j++)
               {
                    ELEM(augmented_matrix, k, j) -= factor * ELEM(augmented_matrix, i, j);
               }
          }
     }
     // oddelenie inverznej casti
     for (i = 0; i < dimension; i++)
     {
          for (j = 0; j < dimension; j++)
          {
               ELEM(inverse_matrix, i, j) = ELEM(augmented_matrix, i, j + dimension);
          }
     }
     mat_destroy(augmented_matrix);
     return inverse_matrix;
}

char mat_division(MAT *a, MAT *b, MAT *c)
{
     MAT *inverse_a;
     unsigned int i, j, k;

     if (b->cols != a->rows)
     {
          return FAILURE;
     }

     inverse_a = mat_invert(a);
     if (inverse_a == NULL)
     {
          return FAILURE;
     }

     mat_zero(c);

     for (i = 0; i < b->rows; i++)
     {

          for (j = 0; j < inverse_a->cols; j++)
          {

               for (k = 0; k < b->cols; k++)
               {
                    ELEM(c, i, j) += ELEM(b, i, k) * ELEM(inverse_a, k, j);
               }
          }
     }
     mat_destroy(inverse_a);
     return SUCCESS;
}

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
          ELEM(coordinatesS, i, 0) = R * cos(Brad) * cos(Lrad); // ELEM je element matice na prislusnom mieste
          ELEM(coordinatesS, i, 1) = R * cos(Brad) * sin(Lrad);
          ELEM(coordinatesS, i, 2) = R * sin(Brad);

          // X koordinaty
          ELEM(coordinatesX, i, 0) = (R+alt) * cos(Brad) * cos(Lrad);
          ELEM(coordinatesX, i, 1) = (R+alt) * cos(Brad) * sin(Lrad);
          ELEM(coordinatesX, i, 2) = (R+alt) * sin(Brad);

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

void calculate_Aij(MAT *A, MAT *coordinatesX, MAT *coordinatesS, MAT *coordinatesE, int n)
{
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               double p1[3] = {ELEM(coordinatesX, i, 0), ELEM(coordinatesX, i, 1), ELEM(coordinatesX, i, 2)};
               double p2[3] = {ELEM(coordinatesS, j, 0), ELEM(coordinatesS, j, 1), ELEM(coordinatesS, j, 2)};
               double dx,dy,dz, rij;
               dx= p1[0] - p2[0];
               dy= p1[1] - p2[1];
               dz= p1[2] - p2[2];
               rij= sqrt(dx*dx + dy*dy + dz*dz);
               double dotProduct;
               dotProduct = dx * ELEM(coordinatesE, i, 0) + dy * ELEM(coordinatesE, i, 1) + dz * ELEM(coordinatesE, i, 2);
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
char mat_multiply(MAT *a, MAT *b, MAT *c)
{
     if (a->cols != b->rows)
     {
          return FAILURE; // Ensure valid multiplication
     }

     // inicializácia matice na výsup
     mat_zero(c); // vynulovanie c

     for (unsigned int i = 0; i < a->rows; i++)
     {
          for (unsigned int j = 0; j < b->cols; j++)
          {
               for (unsigned int k = 0; k < a->cols; k++)
               {
                    ELEM(c, i, j) += ELEM(a, i, k) * ELEM(b, k, j);
               }
          }
     }
     return SUCCESS; // operácia prebehla úspešne
}

void compute_distance_vectors(MAT *coordinatesX, MAT *coordinatesS, MAT *distanceVectors, int n)
{
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               for (int k = 0; k < 3; k++) // Assuming 3D coordinates (x, y, z)
               {
                    // Calculate the difference between corresponding coordinates
                    ELEM(distanceVectors, i, j * 3 + k) = ELEM(coordinatesX, i, k) - ELEM(coordinatesS, j, k);
               }
          }
     }
}

// funkcia na nasobenie matice vektorom, helper function na solver
void mat_vec_mult(MAT *A, MAT *x, MAT *result) {
     mat_zero(result); // Initialize the result vector with zero values
     for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->cols; j++) {
               ELEM(result, i, 0) += ELEM(A, i, j) * ELEM(x, j, 0);
          }
     }
}

// helper function na solver, odpocitanie 2 vektorov
void vec_subtract(MAT *v1, MAT *v2, MAT *result) {
     for (int i = 0; i < v1->rows; i++) {
          ELEM(result, i, 0) = ELEM(v1, i, 0) - ELEM(v2, i, 0);
     }
}

// dot product 2 vektorov, helper function na solver
double vec_dot_product(MAT *v1, MAT *v2) {
     double dot_product = 0.0;
     for (int i = 0; i < v1->rows; i++) {
           dot_product += ELEM(v1, i, 0) * ELEM(v2, i, 0);
     }
     return dot_product;
}

//norma vektora, helper function na solver
double vec_norm(MAT *v) {
     return sqrt(vec_dot_product(v, v));
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

     load_data("/Users/ninalackovicova/Downloads/BL-3602.dat", B, L, H, dg, f, n); // Úprava názvu súboru podľa potreby

     double B_vals[n], L_vals[n];
     for (int i = 0; i < n; i++)
     {
          B_vals[i] = ELEM(B, i, 0); // Get B values
          L_vals[i] = ELEM(L, i, 0); // Get L values
     }

     // funguju funckie multiply a division? TEST
     // K{{1, 3}, {5, 7}} . V{{1, 2}, {6, 7}}  = P{{19, 23}, {47, 59}}
     // MULTIPLY
     /* MAT *K = mat_create_with_type(2, 2);
     ELEM(K, 0, 0) = 1;
     ELEM(K, 0, 1) = 3;
     ELEM(K, 1, 0) = 5;
     ELEM(K, 1, 1) = 7;
     MAT *V = mat_create_with_type(2, 2);
     ELEM(V, 0, 0) = 1;
     ELEM(V, 0, 1) = 2;
     ELEM(V, 1, 0) = 6;
     ELEM(V, 1, 1) = 7;
     MAT *P = mat_create_with_type(2, 2);
     mat_multiply(K, V, P);
     mat_print(P); */

     // Výpočet súradníc
     calculate_coordinates(coordinatesS, coordinatesX, coordinatesE, B_vals, L_vals, n, alt);

     // Výpočet matice vzdialeností
     //calculate_distance_matrix(distanceMatrix, coordinatesX, coordinatesS, n);

     // distanceVectors
     //compute_distance_vectors(coordinatesX, coordinatesS, distanceVectors, n);
     // Výpočet matice Aij

     calculate_Aij(A,  coordinatesX, coordinatesS,coordinatesE, n);

     // Riešenie pre alpha
     MAT *alpha = mat_create_with_type(n, 1);
     MAT *dGMarray = mat_create_with_type(n, 1);

     // naplnenie dGMarray
     for (int i = 0; i < n; i++)
     {
          // dGMarray ma n rows a 1 column
          // ELEM(rows, columns)
          ELEM(dGMarray, i, 1) = dGM;
     }

     // výpočet Gij (musí mať rozmer 1xn), lebo alpha má nx1
     MAT *Gij = mat_create_with_type(1, n);
     // VÝPOČET Gij = 1/(4Pi*rij)
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               // dGMarray ma n rows a 1 column
               // ELEM(rows, columns)
               double rij = ELEM(distanceMatrix, i, j);
               ELEM(Gij, 1, i) = 1 / ((4 * M_PI) * rij);
          }
     }

     // keď neprejde mat_division, kvôli rozmerom, môže invertovať matticu A a vynasobiť ju nasledovne alpha = A(^-1)*dGMarray //DRUHÉ MOŽNÉ RIEŠENIE
     /* if (mat_invert(A) == SUCCESS)
     {
          printf("Alpha vector:\n");
          // mat_print(alpha); //
     }
     else
     {
          printf("Nepodarilo sa vypočítať alphu.\n");
     } */

     // Linerany system A * alpha = dgM, alpha = A^(-1) * dGMarray, A sa v solv. invertuje
     if (mat_division(A, dGMarray, alpha) == SUCCESS)
     {
          printf("Alpha vector:\n");
          // mat_print(alpha); //
     }
     else
     {
          printf("Nepodarilo sa vypočítať alphu.\n");
     }

     // vypísanie dGMarray

     printf("Prvých 10 prvkov Aij:\n");
     for (int i = 0; i < 5; i++)
     {
          for (int j = 0; j < 5; j++)
          {
               printf("%.10f\n", ELEM(A, i, j));

          }
     }

     printf("Prvých 10 prvkov dGMarray:\n");
     for (int i = 0; i < 10; i++)
     {
          printf("%.10f\n", ELEM(dGMarray, i, 1));
     }

     // Free the allocated memory

     // Calculate the final vector u (this assumes some operation with alpha)
     // MAT *u = mat_create_with_type(n, 1);
     // Assuming u = A * alpha for demonstration (replace with your actual calculation)
     // mat_multiply(u, A, alpha); // If you want to compute u = A * alpha

     // Print the result
     // printf("Final vector u:\n");
     // mat_print(u);

     MAT *u1 = mat_create_with_type(1, 1);
     // mat_division(dGMarray, alpha, u1);
     // Násobenie VÝSLEDOK

     mat_multiply(Gij, alpha, u1);
     printf("Final vector u:\n");

     mat_print(u1);

     // uvolenie alokovanej pamäte
     mat_destroy(coordinatesS);
     mat_destroy(coordinatesX);
     mat_destroy(coordinatesE);
     mat_destroy(distanceMatrix);
     mat_destroy(distanceVectors);
     mat_destroy(A);
     mat_destroy(alpha);
     mat_destroy(dGMarray);
     mat_destroy(u1);
     mat_destroy(B);
     mat_destroy(L);
     mat_destroy(H);
     mat_destroy(dg);
     mat_destroy(f);

     return 0;
}
