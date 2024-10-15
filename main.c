
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

double Rtotal(double altitude)
{
     return (R + alt);
}

void calculate_coordinates(MAT *coordinatesS, MAT *coordinatesX, MAT *coordinatesE, double *B, double *L, int n, double altitude)
{
     for (int i = 0; i < n; i++)
     {
          double Brad = rad(B[i]);
          double Lrad = rad(L[i]);

          // S coordinates
          ELEM(coordinatesS, i, 0) = Rtotal(0) * cos(Brad) * cos(Lrad);
          ELEM(coordinatesS, i, 1) = Rtotal(0) * cos(Brad) * sin(Lrad);
          ELEM(coordinatesS, i, 2) = Rtotal(0) * sin(Brad);

          // X coordinates
          ELEM(coordinatesX, i, 0) = Rtotal(altitude) * cos(Brad) * cos(Lrad);
          ELEM(coordinatesX, i, 1) = Rtotal(altitude) * cos(Brad) * sin(Lrad);
          ELEM(coordinatesX, i, 2) = Rtotal(altitude) * sin(Brad);

          // E coordinates
          ELEM(coordinatesE, i, 0) = cos(Brad) * cos(Lrad);
          ELEM(coordinatesE, i, 1) = cos(Brad) * sin(Lrad);
          ELEM(coordinatesE, i, 2) = sin(Brad);
     }
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