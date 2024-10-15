
#include <stdio.h>
#include "mat.h"

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