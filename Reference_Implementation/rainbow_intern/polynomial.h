unsigned i4_choose(unsigned n, unsigned k);

unsigned i4_fall(unsigned x, unsigned n);

unsigned i4_max(unsigned i1, unsigned i2);

unsigned i4_min(unsigned i1, unsigned i2);

void i4vec_concatenate(unsigned n1, const unsigned a[], unsigned n2, const unsigned b[], unsigned c[]);

void i4vec_permute(unsigned n, unsigned p[], unsigned a[]);

unsigned *i4vec_sort_heap_index_a(unsigned n, unsigned a[]);

unsigned i4vec_sum(unsigned n, unsigned a[]);

void mono_next_grlex(unsigned m, unsigned x[]);

unsigned mono_rank_grlex(unsigned m, unsigned x[]);

void mono_total_next_grlex(unsigned m, unsigned n, unsigned x[]);

unsigned *mono_unrank_grlex(unsigned m, unsigned rank);

unsigned mono_upto_enum(unsigned m, unsigned n);

double *mono_value(unsigned m, unsigned n, unsigned f[], double x[]);

void perm_check0(unsigned n, unsigned p[]);

void perm_check1(unsigned n, unsigned p[]);

void polynomial_add(unsigned o1, const unsigned char c1[], const unsigned e1[], unsigned o2, const unsigned char c2[],
                    const unsigned e2[], unsigned *o, unsigned char c[], unsigned offset, unsigned e[]);

void polynomial_axpy(double s, unsigned o1, double c1[], unsigned e1[], unsigned o2, double c2[],
                     unsigned e2[], unsigned *o, double c[], unsigned offset, unsigned e[]);

void polynomial_compress(unsigned o1, unsigned char c1[], unsigned c1_offset, unsigned e1[], unsigned *o2,
                         unsigned char c2[],
                         unsigned c2_offset,
                         unsigned e2[]);

void polynomial_dif(unsigned m, unsigned o1, double c1[], unsigned e1[], unsigned dif[],
                    unsigned *o2, double c2[], unsigned e2[]);

void polynomial_mul(unsigned o1, const unsigned char c1[], unsigned c1_offset, unsigned e1[], unsigned o2,
                    const unsigned char c2[],
                    unsigned c2_offset,
                    unsigned e2[], unsigned *o, unsigned char c[], unsigned c_offset, unsigned e[]);

void polynomial_prunsigned(unsigned o, const unsigned char c[], unsigned gf16_offset, unsigned e[], char *title);

void polynomial_scale(double s, unsigned m, unsigned o1, double c1[], unsigned e1[]);

void polynomial_sort(unsigned o, unsigned char c[], unsigned offset, unsigned e[]);

double *polynomial_value(unsigned m, unsigned o, double c[], unsigned e[], unsigned nx,
                         double x[]);

void
r8vec_concatenate(unsigned n1, unsigned char a[], unsigned n2, unsigned char b[], unsigned char c[], unsigned offset);

void r8vec_permute(unsigned n, unsigned p[], unsigned char a[], unsigned offset);