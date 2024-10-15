

#define ELEM(mat, i, j) ((mat)->elem[(i) * (mat)->cols + j])
#define SUCCESS 1
#define FAILURE 0

typedef struct
{
     unsigned int rows;
     unsigned int cols;
     float *elem;
} MAT;

MAT *mat_create_with_type(unsigned int rows, unsigned int cols);
MAT *mat_zero(MAT *null);
MAT *mat_create_by_file(char *filename);
char mat_save(MAT *mat, char *filename);
void mat_destroy(MAT *mat);
void mat_unit(MAT *mat);
void mat_random(MAT *mat);
void mat_print(MAT *mat);
MAT *mat_invert(MAT *input_matrix);
char mat_division(MAT *a, MAT *b, MAT *c);
