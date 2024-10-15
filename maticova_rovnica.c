#include <stdio.h>
#include <stdlib.h>
#define ELEM(mat, i, j) ((mat)->elem[(i) * (mat)->cols + j])
#define SUCCESS 1
#define FAILURE 0

typedef struct
{
     unsigned int rows;
     unsigned int cols;
     float *elem;
} MAT;

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
               printf("%.2f\t", ELEM(mat, i, j));
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

void main()
{
     MAT *a, *b, *m;
     MAT *nacitana;

     a = mat_create_with_type(2, 2);
     b = mat_create_with_type(2, 2);
     m = mat_create_with_type(2, 2);

     ELEM(a, 0, 0) = 1;
     ELEM(a, 0, 1) = 2;
     ELEM(a, 1, 0) = 3;
     ELEM(a, 1, 1) = 4;

     ELEM(b, 0, 0) = 5;
     ELEM(b, 0, 1) = 3;
     ELEM(b, 1, 0) = 1;
     ELEM(b, 1, 1) = 7;

     printf("A = \n");
     mat_print(a);
     printf("B =\n");
     mat_print(b);

     if (mat_division(a, b, m) == SUCCESS)
     {
          printf("B/A = \n");
          mat_print(m);
     }
     else
     {
          printf("Fail.\n");
     }

     // ulozenie matice
     if (mat_save(a, "mat_a.dat") == SUCCESS)
     {
          printf("Matica bola uspesne ulozena.\n");
     }
     else
     {
          printf("Ulozenie sa nepodarilo.\n");
     }

     // nacitanie matice
     nacitana = mat_create_by_file("mat_a.dat");
     if (nacitana != NULL)
     {
          printf("Matica bola uspesne nacitana.\n");
          mat_print(nacitana);
     }
     else
     {
          printf("Nacitanie zlyhalo.\n");
     }

     mat_destroy(a);
     mat_destroy(b);
     mat_destroy(m);
     mat_destroy(nacitana);
}