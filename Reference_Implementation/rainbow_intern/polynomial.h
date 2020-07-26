unsigned i4_choose(int n, int k);

int i4_max(int i1, int i2);

int i4_min(int i1, int i2);

void i4vec_concatenate(unsigned n1, const unsigned a[], unsigned n2, const unsigned b[], unsigned c[]);

void i4vec_permute(unsigned n, unsigned p[], unsigned a[]);

unsigned *i4vec_sort_heap_index_a(unsigned n, const unsigned int *a);

int i4vec_sum(unsigned n, const unsigned int *a);

unsigned mono_rank_grlex(int m, unsigned *x);

void mono_unrank_grlex(int m, unsigned rank, unsigned *dest_f);

unsigned char mono_value(unsigned f[], unsigned char x[]);

void perm_check0(unsigned n, const unsigned int *p);

void
polynomial_add(unsigned char *destSummand, unsigned dest_offset, unsigned dest_grade, const unsigned char *summand,
               unsigned summand_offset, unsigned summand_o, const unsigned int *summand_e);

void polynomial_compress(unsigned o1, unsigned char c1[], const unsigned int e1[], unsigned *o2, unsigned char c2[],
                         unsigned int e2[]);

void
polynomial_mul(const unsigned char *factor_A, unsigned A_offset, unsigned A_grade, const unsigned char *factor_B,
               unsigned B_offset, unsigned B_grade, unsigned char *C, unsigned *C_o, unsigned int C_e[]);

void polynomial_print(unsigned o, const unsigned char *c, unsigned gf16_offset, const unsigned *e, char *title);

void polynomial_sort(unsigned o, unsigned char c[], unsigned int e[]);

unsigned char
polynomial_value(unsigned o, const unsigned char *c, unsigned offset, unsigned const e[], unsigned char *x);

void
r8vec_concatenate(unsigned n1, const unsigned char a[], unsigned n2, const unsigned char b[], unsigned char c[],
                  unsigned offset);

void r8vec_permute(unsigned n1, unsigned int p1[], unsigned char a[]);