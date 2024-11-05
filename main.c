
#include <stdio.h>
#include <math.h>
#include "mat.h"
#include <stdlib.h>

#define R 6378
#define N 8102
#define GM 398600.5
#define alt 230
#define dGM 0.0000270141
#define TOL 1e-9 // Tolerancia na konvergenciu

// Funckia na násobenie matice a k nej transponovanej matice
// nevytváram explicitne AT

MAT *mat_multiply_with_transpose(MAT *a)
{
     // musia sedieť dimenzie, resp. najprv, či nie sú nulové
     if (!a || a->rows == 0 || a->cols == 0)
     {
          // vysledna matica bude vždy štvorcová
          return NULL; // matica nie je dobre zadaná
     }

     // alokovanie matice
     MAT *result = mat_create_with_type(a->rows, a->rows);
     if (!result)
     {
          return NULL; // zle alokovane
     }

     //  A * AT
     // vysledna[i][j] = sum (od k=0, po n-1) a[i][k]xa[j][k]
     for (unsigned int i = 0; i < a->rows; i++)
     {
          for (unsigned int j = 0; j < a->rows; j++)
          {
               // vynulovanie
               ELEM(result, i, j) = 0.0f;
               // dot product
               for (unsigned int k = 0; k < a->cols; k++)
               {
                    ELEM(result, i, j) += ELEM(a, i, k) * ELEM(a, j, k);
               }
          }
     }

     return result;
}

MAT *mat_create_with_type(unsigned int rows, unsigned int cols)
{
     MAT *mat = (MAT *)malloc(sizeof(MAT));
     if (!mat)
     {
          fprintf(stderr, "Alokacia pamate pre maticu zlyhala, strukturalne.\n");
          return NULL;
     }

     if (mat == NULL)
     {
          return NULL;
     }
     mat->rows = rows;
     mat->cols = cols;
     mat->elem = (float *)malloc(rows * cols * sizeof(float));
     if (mat->elem == NULL)
     {
          fprintf(stderr, "Alokacia pamate pre maticu zlyhala, nepodarilo sa alokovat elementy.\n");
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
          ELEM(coordinatesX, i, 0) = (R + alt) * cos(Brad) * cos(Lrad);
          ELEM(coordinatesX, i, 1) = (R + alt) * cos(Brad) * sin(Lrad);
          ELEM(coordinatesX, i, 2) = (R + alt) * sin(Brad);

          // E koordinaty
          ELEM(coordinatesE, i, 0) = cos(Brad) * cos(Lrad);
          ELEM(coordinatesE, i, 1) = cos(Brad) * sin(Lrad);
          ELEM(coordinatesE, i, 2) = sin(Brad);
     }
}

/*double distance(double *x1, double *x2)
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
}*/

// pocitanie vzdialenosti
void calculate_distance_matrix(MAT *distanceMatrix, MAT *coordinatesX, MAT *coordinatesS, int n)
{
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               double dx = ELEM(coordinatesX, i, 0) - ELEM(coordinatesS, j, 0);
               double dy = ELEM(coordinatesX, i, 1) - ELEM(coordinatesS, j, 1);
               double dz = ELEM(coordinatesX, i, 2) - ELEM(coordinatesS, j, 2);
               ELEM(distanceMatrix, i, j) = sqrt(dx * dx + dy * dy + dz * dz);
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
               double dx, dy, dz, rij;
               dx = p1[0] - p2[0];
               dy = p1[1] - p2[1];
               dz = p1[2] - p2[2];
               rij = sqrt(dx * dx + dy * dy + dz * dz);
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
void mat_vec_mult(MAT *A, MAT *x, MAT *result)
{

     if (A->cols != x->rows)
     {
          fprintf(stderr, "Error: zle dimenzie na funckiu Matica krat vektor\n");
          return;
     }
     mat_zero(result); // Initialize the result vector with zero values
     for (int i = 0; i < A->rows; i++)
     {
          for (int j = 0; j < A->cols; j++)
          {
               ELEM(result, i, 0) += ELEM(A, i, j) * ELEM(x, j, 0);
          }
     }
}

// helper function na solver, odpocitanie 2 vektorov
void vec_subtract(MAT *v1, MAT *v2, MAT *result)
{
     if (v1->rows != v2->rows)
     {
          fprintf(stderr, "Error: vektory nemaju rovnaku dlzku\n");
          return;
     }
     for (int i = 0; i < v1->rows; i++)
     {
          ELEM(result, i, 0) = ELEM(v1, i, 0) - ELEM(v2, i, 0);
     }
}

// dot product 2 vektorov, helper function na solver
double vec_dot_product(MAT *v1, MAT *v2)
{
     if (v1->rows != v2->rows)
     {
          fprintf(stderr, "Error: vektory musia mat rovnaku dlzku\n");
          return 0.0; // or handle the error appropriately
     }
     double dot_product = 0.0;
     for (int i = 0; i < v1->rows; i++)
     {
          dot_product += ELEM(v1, i, 0) * ELEM(v2, i, 0);
     }
     return dot_product;
}

// norma vektora, helper function na solver
double vec_norm(MAT *v)
{
     return sqrt(vec_dot_product(v, v));
}

int BiCGSTAB(MAT *A, MAT *b, MAT *x)
{
     int n = A->rows;
     // vytvorenie pomocných vektorov
     MAT *r = mat_create_with_type(n, 1);
     MAT *r_hat = mat_create_with_type(n, 1);
     MAT *v = mat_create_with_type(n, 1);
     MAT *p = mat_create_with_type(n, 1);
     MAT *s = mat_create_with_type(n, 1);
     MAT *t = mat_create_with_type(n, 1);
     MAT *Ap = mat_create_with_type(n, 1);

     // Skalárne premenné pre BiCGSTAB algoritmus
     double rho = 0.0, alpha = 1.0, omega = 1.0;
     double beta;
     double rho_prev = rho; // nastavenie pred 1.iteraciou

     // Výpočet počiatočného rezidua r(0) = b - A * x(0)
     mat_vec_mult(A, x, Ap); // Ap = A * x
     vec_subtract(b, Ap, r); // r = b - Ap

     // pociatocne reziduum
     /* printf("Pociatocne reziduum r:\n");
     for (int i = 0; i < n; i++)
     {
          printf("%f ", ELEM(r, i, 0));
     }
     printf("\n"); */

     // Nastavenie r_hat = r (predpokladáme, že r_hat je kópia počiatočného rezidua)
     mat_zero(r_hat);
     for (int i = 0; i < n; i++)
     {
          ELEM(r_hat, i, 0) = ELEM(r, i, 0); // r_hat = r
     }

     // Debugging pre r_hat
     printf("Pociatocne r_hat:\n");
     for (int i = 0; i < n; i++)
     {
          printf("%f ", ELEM(r_hat, i, 0));
     }
     printf("\n");

     // BiCGSTAB iterations
     for (int iter = 0; iter < n; iter++)
     {
          // Výpočet ρ(i-1) = r_hat^T * r(i-1)
          rho_prev = rho;
          rho = vec_dot_product(r_hat, r); //// ρ(i-1) = r_hat^T * r

          printf("Iteration: %d, rho: %.10f\n", iter, rho);

          if (fabs(rho) < TOL)
          {
               printf("Zlyhanie v BiCGSTAB: rho je príliš malé\n");
               break;
          }

          // Ak i = 1, nastav p(1) = r(0)
          if (iter == 0)
          {
               for (int i = 0; i < n; i++)
               {
                    ELEM(p, i, 0) = ELEM(r, i, 0); // p0 = r0
               }
               // inak smerový vektor p(i) = r(i-1) + β(i-1)(p(i-1) - ω(i-1)v(i-1))
          }
          else
          {
               beta = (rho / rho_prev) * (alpha / omega);
               for (int i = 0; i < n; i++)
               {
                    ELEM(p, i, 0) = ELEM(r, i, 0) + beta * (ELEM(p, i, 0) - omega * ELEM(v, i, 0));
               }
          }

          // v = A * p
          mat_vec_mult(A, p, v);

          // α(i) = ρ(i-1) / (r_hat^T * v(i))
          alpha = rho / vec_dot_product(r_hat, v);

          // s = r(i-1) - α(i) * v(i)
          for (int i = 0; i < n; i++)
          {
               ELEM(s, i, 0) = ELEM(r, i, 0) - alpha * ELEM(v, i, 0);
          }

          // Skontroluj normu vektora s; ak je dostatočne malá, ukonči iteráciu
          if (vec_norm(s) < TOL)
          {
               for (int i = 0; i < n; i++)
               {
                    ELEM(x, i, 0) += alpha * ELEM(p, i, 0); // x(i) = x(i-1) + α(i) * p(i)
               }
               break;
          }

          // t = A * s
          mat_vec_mult(A, s, t);

          // ω(i) = (t^T * s) / (t^T * t)
          omega = vec_dot_product(t, s) / vec_dot_product(t, t);

          // Aktualizácia riešenia x(i) = x(i-1) + α(i) * p(i) + ω(i) * s(i)
          for (int i = 0; i < n; i++)
          {
               ELEM(x, i, 0) += alpha * ELEM(p, i, 0) + omega * ELEM(s, i, 0);
          }

          // Update r(i) = s - ω(i) * t
          for (int i = 0; i < n; i++)
          {
               ELEM(r, i, 0) = ELEM(s, i, 0) - omega * ELEM(t, i, 0);
          }

          // Checkni normu nového rezidua r; ak je dostatočne malé, ukonči iteráciu
          if (vec_norm(r) < TOL)
          {
               break;
          }

          // ak ω(i) = 0 -> delenie nulov -> zlyhanie
          if (fabs(omega) < TOL)
          {
               printf("Breakdown in BiCGSTAB: omega is too small\n");
               break;
          }
     }

     // Clean up
     mat_destroy(r);
     mat_destroy(r_hat);
     mat_destroy(v);
     mat_destroy(p);
     mat_destroy(s);
     mat_destroy(t);
     mat_destroy(Ap);

     return 0;
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

     load_data("/Users/hannah/Desktop/GOCE-2/BL-8102.dat", B, L, H, dg, f, n); // Úprava názvu súboru podľa potreby

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
     // calculate_distance_matrix(distanceMatrix, coordinatesX, coordinatesS, n);

     // distanceVectors
     // compute_distance_vectors(coordinatesX, coordinatesS, distanceVectors, n);
     // Výpočet matice Aij

     calculate_Aij(A, coordinatesX, coordinatesS, coordinatesE, n);

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

     // BiCGSTAB solver
     if (BiCGSTAB(A, dGMarray, alpha) != SUCCESS)
     {
          fprintf(stderr, "BiCGSTAB sa nepodarilo.\n");
          // Handle the failure case
     }

     printf("Prvých 10 prvkov Alphy:\n");
     for (int i = 0; i < 10; i++)
     {
          printf("%.10f\n", ELEM(alpha, i, 1));
     }

     // výpočet Gij (musí mať rozmer 1xn), lebo alpha má nx1
     /* MAT *Gij = mat_create_with_type(1, n);
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
     }  */

     // keď neprejde mat_division, kvôli rozmerom, môže invertovať matticu A a vynasobiť ju nasledovne alpha = A(^-1)*dGMarray //DRUHÉ MOŽNÉ RIEŠENIE
     /* if (mat_invert(A) == SUCCESS)
     {
          printf("Alpha vector:\n");
          // mat_print(alpha); //
     }
     else
     {
          printf("Nepodarilo sa vypočítať alphu.\n");
     }  */

     // Linerany system A * alpha = dgM, alpha = A^(-1) * dGMarray, A sa v solv. invertuje
     /* if (mat_division(A, dGMarray, alpha) == SUCCESS)
     {
          printf("Alpha vector:\n");
          // mat_print(alpha); //
     }
     else
     {
          printf("Nepodarilo sa vypočítať alphu.\n");
     }
 */
     // vypísanie dGMarray

     printf("Prvých 10 prvkov Aij:\n");
     for (int i = 0; i < 5; i++)
     {
          for (int j = 0; j < 5; j++)
          {
               printf("%.10f\n", ELEM(A, i, j));
          }
     }

     /*  printf("Prvých 10 prvkov dGMarray:\n");
      for (int i = 0; i < 10; i++)
      {
           printf("%.10f\n", ELEM(dGMarray, i, 1));
      } */

     // Free the allocated memory

     // Calculate the final vector u (this assumes some operation with alpha)
     // MAT *u = mat_create_with_type(n, 1);
     // Assuming u = A * alpha for demonstration (replace with your actual calculation)
     // mat_multiply(u, A, alpha); // If you want to compute u = A * alpha

     // Print the result
     // printf("Final vector u:\n");
     // mat_print(u);

     /* MAT *u1 = mat_create_with_type(1, 1);
     // mat_division(dGMarray, alpha, u1);
     // Násobenie VÝSLEDOK

     mat_multiply(Gij, alpha, u1);
     printf("Final vector u:\n");

     mat_print(u1); */

     // EFEKTÍVNE NÁSOBENIE MATICE A K NEJ TRASPONOVANEJ MATICE

     // skuska
     MAT *Q = mat_create_with_type(2, 3);
     if (Q)
     {
          // naplnenie matice
          ELEM(Q, 0, 0) = 1.0f;
          ELEM(Q, 0, 1) = 2.0f;
          ELEM(Q, 0, 2) = 3.0f;
          ELEM(Q, 1, 0) = 4.0f;
          ELEM(Q, 1, 1) = 5.0f;
          ELEM(Q, 1, 2) = 6.0f;

          printf("Matica A:\n");
          mat_print(Q);

          // A krat A transponovane
          MAT *result = mat_multiply_with_transpose(Q);
          if (result)
          {
               printf("Výsledok A * A^T:\n");
               mat_print(result);

               // uvolnenie pamate
               free(result->elem);
               free(result);
          }

          free(Q->elem);
          free(Q);
     }

     // uvolenie alokovanej pamäte
     mat_destroy(coordinatesS);
     mat_destroy(coordinatesX);
     mat_destroy(coordinatesE);
     mat_destroy(distanceMatrix);
     mat_destroy(distanceVectors);
     mat_destroy(A);
     mat_destroy(alpha);
     mat_destroy(dGMarray);
     // mat_destroy(u1);
     mat_destroy(B);
     mat_destroy(L);
     mat_destroy(H);
     mat_destroy(dg);
     mat_destroy(f);

     return 0;
}
