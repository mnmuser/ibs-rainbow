unsigned i4_choose(unsigned n, unsigned k);

unsigned i4_max(unsigned i1, unsigned i2);

unsigned i4_min(unsigned i1, unsigned i2);

void i4vec_concatenate(unsigned n1, const unsigned a[], unsigned n2, const unsigned b[], unsigned c[]);

void i4vec_permute(unsigned n, unsigned p[], unsigned a[]);

unsigned *i4vec_sort_heap_index_a(unsigned n, unsigned a[]);

unsigned i4vec_sum(unsigned n, unsigned a[]);

void mono_next_grlex(unsigned m, unsigned x[]);

unsigned mono_rank_grlex(unsigned m, unsigned x[]);

unsigned *mono_unrank_grlex(unsigned m, unsigned rank);

unsigned char mono_value(unsigned f[], unsigned char x[]);

void perm_check0(unsigned n, unsigned p[]);

void polynomial_add(unsigned char *dest, unsigned dest_offset, unsigned *dest_o, unsigned int dest_e[], unsigned A_o,
                    unsigned char *summand_A, unsigned A_offset, const unsigned int A_e[], unsigned B_o,
                    const unsigned char *summand_B, unsigned B_offset, const unsigned int B_e[]);

void polynomial_compress(unsigned o1, unsigned char c1[], unsigned c1_offset, unsigned e1[], unsigned *o2,
                         unsigned char c2[],
                         unsigned c2_offset,
                         unsigned e2[]);

void polynomial_mul(unsigned o1, const unsigned char c1[], unsigned c1_offset, const unsigned e1[], unsigned o2,
                    const unsigned char c2[],
                    unsigned c2_offset,
                    const unsigned e2[], unsigned *o, unsigned char c[], unsigned c_offset, unsigned e[]);

void polynomial_print(unsigned o, const unsigned char *c, unsigned gf16_offset, const unsigned *e, char *title);

void polynomial_sort(unsigned o, unsigned char c[], unsigned offset, unsigned e[]);

unsigned char
polynomial_value(unsigned o, const unsigned char *c, unsigned offset, unsigned const e[], unsigned char *x);

void
r8vec_concatenate(unsigned n1, const unsigned char a[], unsigned n2, const unsigned char b[], unsigned char c[],
                  unsigned offset);

void r8vec_permute(unsigned n, unsigned p[], unsigned char a[], unsigned offset);