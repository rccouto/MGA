char *Char_alloc(char *name2);
char *Char_free(char *name);
int *Int_alloca (int n);
int *Int_freea (int *v);
double *Real_alloca (int n);
double *Real_freea (double *v);
int **Int_allocm (int m, int n);
int **Int_freem (int m, int n, int **v);
double **Real_allocm (int m, int n);
double **Real_freem (int m, int n, double **v);