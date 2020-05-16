int i4_choose(int n, int k);

int i4_fall(int x, int n);

int i4_max(int i1, int i2);

int i4_min(int i1, int i2);

void i4vec_concatenate(int n1, int a[], int n2, int b[], int c[]);

void i4vec_permute(int n, int p[], int a[]);

int *i4vec_sort_heap_index_a(int n, int a[]);

int i4vec_sum(int n, int a[]);

void mono_next_grlex(int m, int x[]);

int mono_rank_grlex(int m, int x[]);

void mono_total_next_grlex(int m, int n, int x[]);

int *mono_unrank_grlex(int m, int rank);

int mono_upto_enum(int m, int n);

double *mono_value(int m, int n, int f[], double x[]);

void perm_check0(int n, int p[]);

void perm_check1(int n, int p[]);

void polynomial_add(int o1, unsigned char c1[], int e1[], int o2, unsigned char c2[],
                    int e2[], int *o, unsigned char c[], unsigned offset, int e[]);

void polynomial_axpy(double s, int o1, double c1[], int e1[], int o2, double c2[],
                     int e2[], int *o, double c[], unsigned offset, int e[]);

void polynomial_compress(int o1, unsigned char c1[], int e1[], int *o2, unsigned char c2[],
                         int e2[]);

void polynomial_dif(int m, int o1, double c1[], int e1[], int dif[],
                    int *o2, double c2[], int e2[]);

void polynomial_mul(int o1, const unsigned char c1[], unsigned c1_offset, int e1[], int o2, const unsigned char c2[],
                    unsigned c2_offset,
                    int e2[], int *o, unsigned char c[], int e[]);

void polynomial_print(int m, int o, const unsigned char c[], unsigned gf16_offset, int e[], char *title);

void polynomial_scale(double s, int m, int o1, double c1[], int e1[]);

void polynomial_sort(int o, unsigned char c[], int e[]);

double *polynomial_value(int m, int o, double c[], int e[], int nx,
                         double x[]);

void r8vec_concatenate(int n1, unsigned char a[], int n2, unsigned char b[], unsigned char c[], unsigned offset);

void r8vec_permute(int n, int p[], unsigned char a[]);