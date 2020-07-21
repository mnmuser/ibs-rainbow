#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <memory.h>

#include "polynomial.h"
#include "rainbow_blas.h"
#include "rainbow_keypair_computation.h"


/*
 * Code from John Burkhardt, GNU LGPL licensed.
 * https://people.sc.fsu.edu/~jburkardt/c_src/polynomial/polynomial.html
 *
 * Needed for Operations on the quartic id-polynomials.
 * Adapted for GF16 (and GF256).
 *
 * Modified in May 2020.
 */

/******************************************************************************/

unsigned i4_choose(int n, int k)

/******************************************************************************/
/*
  Purpose:

    I4_CHOOSE computes the binomial coefficient C(N,K).

  Discussion:

    The value is calculated in such a way as to avoid overflow and
    roundoff.  The calculation is done in integer arithmetic.

    The formula used is:

      C(N,K) = N! / ( K! * (N-K)! )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 May 2008

  Author:

    John Burkardt

  Reference:

    ML Wolfson, HV Wright,
    Algorithm 160:
    Combinatorial of M Things Taken N at a Time,
    Communications of the ACM,
    Volume 6, Number 4, April 1963, page 161.

  Parameters:

    Input, unsigned N, K, are the values of N and K.

    Output, unsigned I4_CHOOSE, the number of combinations of N
    things taken K at a time.
*/
{
    int i;
    int mn;
    int mx;
    unsigned value;

    mn = i4_min(k, n - k);

    if (mn < 0) {
        value = 0;
    } else if (mn == 0) {
        value = 1;
    } else {
        mx = i4_max(k, n - k);
        value = mx + 1;

        for (i = 2; i <= mn; i++) {
            value = (value * (mx + i)) / i;
        }
    }

    return value;
}

/******************************************************************************/

int i4_max(int i1, int i2)

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, unsigned I1, I2, are two ints to be compared.

    Output, unsigned I4_MAX, the larger of I1 and I2.
*/
{
    int value;

    if (i2 < i1) {
        value = i1;
    } else {
        value = i2;
    }
    return value;
}

/******************************************************************************/

int i4_min(int i1, int i2)

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, unsigned I1, I2, two ints to be compared.

    Output, unsigned I4_MIN, the smaller of I1 and I2.
*/
{
    int value;

    if (i1 < i2) {
        value = i1;
    } else {
        value = i2;
    }
    return value;
}

/******************************************************************************/

void i4vec_concatenate(unsigned n1, const unsigned a[], unsigned n2, const unsigned b[], unsigned c[])

/******************************************************************************/
/*
  Purpose:

    I4VEC_CONCATENATE concatenates two I4VEC's.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, unsigned N1, the number of entries in the first vector.

    Input, unsigned A[N1], the first vector.

    Input, unsigned N2, the number of entries in the second vector.

    Input, unsigned B[N2], the second vector.

    Output, unsigned I4VEC_CONCATENATE[N1+N2], the concatenated vector.
*/
{
    unsigned i;

    for (i = 0; i < n1; i++) {
        c[i] = a[i];
    }
    for (i = 0; i < n2; i++) {
        c[n1 + i] = b[i];
    }
}

/******************************************************************************/

void i4vec_permute(unsigned n1, unsigned p1[], unsigned a[])

/******************************************************************************/
/*
  Purpose:

    I4VEC_PERMUTE permutes an I4VEC in place.

  Discussion:

    An I4VEC is a vector of I4's.

    This routine permutes an array of integer "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.

  Example:

    Input:

      N = 5
      P = (   1,   3,   4,   0,   2 )
      A = (   1,   2,   3,   4,   5 )

    Output:

      A    = (   2,   4,   5,   1,   3 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects.

    Input, int P[N], the permutation.  P(I) = J means
    that the I-th element of the output array should be the J-th
    element of the input array.

    Input/output, int A[N], the array to be permuted.
*/
{
    int *p = (int *) p1; //CAST unsigned to signed for this operations
    int n = (int) n1;

    unsigned a_temp;
    int i;
    int iget;
    int iput;
    int istart;

    perm_check0(n, p1);
/*
  In order for the sign negation trick to work, we need to assume that the
  entries of P are strictly positive.  Presumably, the lowest number is 0.
  So temporarily add 1 to each entry to force positivity.
*/
    for (i = 0; i < n; i++) {
        p[i] = p[i] + 1;
    }
/*
  Search for the next element of the permutation that has not been used.
*/
    for (istart = 1; istart <= n; istart++) {
        if (p[istart - 1] < 0) {
            continue;
        } else if (p[istart - 1] == istart) {
            p[istart - 1] = -p[istart - 1];
            continue;
        } else {
            a_temp = a[istart - 1];
            iget = istart;
/*
  Copy the new value into the vacated entry.
*/
            for (;;) {
                iput = iget;
                iget = p[iget - 1];

                p[iput - 1] = -p[iput - 1]; //TODO: hier geht was schief (?)

                if (iget < 1 || n < iget) {
                    fprintf(stderr, "\n");
                    fprintf(stderr, "I4VEC_PERMUTE - Fatal error!\n");
                    fprintf(stderr, "  Entry IPUT = %d of the permutation has\n", iput);
                    fprintf(stderr, "  an illegal value IGET = %d.\n", iget);
                    exit(1);
                }

                if (iget == istart) {
                    a[iput - 1] = a_temp;
                    break;
                }
                a[iput - 1] = a[iget - 1];
            }
        }
    }
/*
  Restore the signs of the entries.
*/
    for (i = 0; i < n; i++) {
        p[i] = -p[i];
    }
/*
  Restore the entries.
*/
    for (i = 0; i < n; i++) {
        p[i] = p[i] - 1;
    }
}

/******************************************************************************/

unsigned *i4vec_sort_heap_index_a(unsigned n, const unsigned *a)

/******************************************************************************/
/*
  Purpose:

    I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(indx(*))

    or explicitly, by the call

      i4vec_permute ( n, indx, a )

    after which a(*) is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, unsigned N, the number of entries in the array.

    Input, unsigned A[N], an array to be index-sorted.

    Output, unsigned I4VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
    unsigned aval;
    unsigned i;
    unsigned *indx;
    unsigned indxt;
    unsigned ir;
    unsigned j;
    unsigned l;

    if (n < 1) {
        return NULL;
    }

    indx = (unsigned *) malloc(n * sizeof(unsigned));

    for (i = 0; i < n; i++) {
        indx[i] = i;
    }

    if (n == 1) {
        indx[0] = indx[0];
        return indx;
    }

    l = n / 2 + 1;
    ir = n;

    for (;;) {

        if (1 < l) {
            l = l - 1;
            indxt = indx[l - 1];
            aval = a[indxt];
        } else {
            indxt = indx[ir - 1];
            aval = a[indxt];
            indx[ir - 1] = indx[0];
            ir = ir - 1;

            if (ir == 1) {
                indx[0] = indxt;
                break;
            }
        }

        i = l;
        j = l + l;

        while (j <= ir) {
            if (j < ir) {
                if (a[indx[j - 1]] < a[indx[j]]) {
                    j = j + 1;
                }
            }

            if (aval < a[indx[j - 1]]) {
                indx[i - 1] = indx[j - 1];
                i = j;
                j = j + j;
            } else {
                j = ir + 1;
            }
        }
        indx[i - 1] = indxt;
    }

    return indx;
}

/******************************************************************************/

int i4vec_sum(unsigned n, const unsigned *a)

/******************************************************************************/
/*
  Purpose:

    I4VEC_SUM sums the entries of an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Example:

    Input:

      A = ( 1, 2, 3, 4 )

    Output:

      I4VEC_SUM = 10

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, unsigned N, the number of entries in the vector.

    Input, unsigned A[N], the vector to be summed.

    Output, unsigned I4VEC_SUM, the sum of the entries of A.
*/
{
    unsigned i;
    int sum;

    sum = 0;
    for (i = 0; i < n; i++) {
        sum = sum + (int) a[i];
    }

    return sum;
}

unsigned mono_rank_grlex(int m, unsigned *x)

/******************************************************************************/
/*
  Purpose:

    MONO_RANK_GRLEX computes the graded lexicographic rank of a monomial.

  Discussion:

    The graded lexicographic ordering is used, over all monomials in 
    M dimensions, for total degree = 0, 1, 2, ...

    For example, if M = 3, the ranking begins:

    Rank  Sum    1  2  3
    ----  ---   -- -- --
       1    0    0  0  0

       2    1    0  0  1
       3    1    0  1  0
       4    1    1  0  1

       5    2    0  0  2
       6    2    0  1  1
       7    2    0  2  0
       8    2    1  0  1
       9    2    1  1  0
      10    2    2  0  0

      11    3    0  0  3
      12    3    0  1  2
      13    3    0  2  1
      14    3    0  3  0
      15    3    1  0  2
      16    3    1  1  1
      17    3    1  2  0
      18    3    2  0  1
      19    3    2  1  0
      20    3    3  0  0

      21    4    0  0  4
      ..   ..   .. .. ..

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, unsigned M, the spatial dimension.
    1 <= M.

    Input, unsigned X[M], the composition.
    For each 1 <= I <= M, we have 0 <= X(I).

    Output, unsigned MONO_RANK_GRLEX, the rank.
*/
{
    int i;
    int j;
    int ks;
    int n;
    int nm;
    int ns;
    unsigned rank;
    int tim1;
    int *xs;
/*
  Ensure that 1 <= M.
*/
    if (m < 1) {
        fprintf(stderr, "\n");
        fprintf(stderr, "MONO_RANK_GRLEX - Fatal error!\n");
        fprintf(stderr, "  M < 1\n");
        exit(1);
    }
/*
  NM = sum ( X )
*/
    nm = i4vec_sum(m, x);
/*
  Convert to KSUBSET format.
*/
    ns = nm + m - 1;
    ks = m - 1;
    xs = (int *) malloc(ks * sizeof(int));
    xs[0] = (int) x[0] + 1;
    for (i = 2; i < m; i++) {
        xs[i - 1] = xs[i - 2] + (int) x[i - 1] + 1;
    }
/*
  Compute the rank.
*/
    rank = 1;

    for (i = 1; i <= ks; i++) {
        if (i == 1) {
            tim1 = 0;
        } else {
            tim1 = xs[i - 2];
        }

        if (tim1 + 1 <= xs[i - 1] - 1) {
            for (j = tim1 + 1; j <= xs[i - 1] - 1; j++) {
                rank = rank + i4_choose(ns - j, ks - i);
            }
        }
    }

    for (n = 0; n < nm; n++) {
        rank = rank + i4_choose(n + m - 1, n);
    }

    free(xs);

    return rank;
}

/******************************************************************************/

void mono_unrank_grlex(int m, unsigned rank, unsigned *dest_f)

/******************************************************************************/
/*
  Purpose:

    MONO_UNRANK_GRLEX computes the monomial of given grlex rank.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2014

  Author:

    John Burkardt

  Parameters:

    Input, unsigned M, the spatial dimension.
    1 <= M.

    Input, unsigned RANK, the rank of the composition.
    1 <= RANK.

    Output, unsigned COMP_UNRANK_GRLEX[M], the composition X of the given rank.
    For each I, 0 <= X[I] <= NM, and 
    sum ( 1 <= I <= M ) X[I] = NM.
*/
{
    int i;
    int j;
    int ks;
    int nm;
    int ns;
    unsigned r;
    unsigned rank1;
    unsigned rank2;

    // removed malloc and free concept (+ 30% performance while debugging)
/*
  Ensure that 1 <= M.
*/
    if (m < 1) {
        fprintf(stderr, "\n");
        fprintf(stderr, "MONO_UNRANK_GRLEX - Fatal error!\n");
        fprintf(stderr, "  M < 1\n");
        exit(1);
    }
/*
  Ensure that 1 <= RANK.
*/
    if (rank < 1) {
        fprintf(stderr, "\n");
        fprintf(stderr, "MONO_UNRANK_GRLEX - Fatal error!\n");
        fprintf(stderr, "  RANK < 1\n");
        exit(1);
    }
/*
  Special case M == 1.
*/
    if (m == 1) {
        dest_f[0] = rank - 1;
        return;
    }
/*
  Determine the appropriate value of NM.
  Do this by adding up the number of compositions of sum 0, 1, 2, 
  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
  gives you the rank of the composition within the set of compositions
  of sum NM.  And that's the number you need in order to do the
  unranking.
*/
    rank1 = 1;
    nm = -1;
    for (;;) {
        nm = nm + 1;
        r = i4_choose(nm + m - 1, nm);
        if (rank < rank1 + r) {
            break;
        }
        rank1 = rank1 + r;
    }

    rank2 = rank - rank1;
/*
  Convert to KSUBSET format.
  Apology: an unranking algorithm was available for KSUBSETS,
  but not immediately for compositions.  One day we will come back
  and simplify all this.
*/
    ks = m - 1;
    ns = nm + m - 1;
    unsigned xs[ks];

    j = 1;

    for (i = 1; i <= ks; i++) {
        r = i4_choose(ns - j, ks - i);

        while (r <= rank2 && 0 < r) {
            rank2 = rank2 - r;
            j = j + 1;
            r = i4_choose(ns - j, ks - i);
        }
        xs[i - 1] = j;
        j = j + 1;
    }
/*
  Convert from KSUBSET format to COMP format.
*/
    dest_f[0] = xs[0] - 1;
    for (i = 2; i < m; i++) {
        dest_f[i - 1] = xs[i - 1] - xs[i - 2] - 1;
    }
    dest_f[m - 1] = ns - xs[ks - 1];
}

/******************************************************************************/

unsigned char mono_value(unsigned f[], unsigned char x[])

/******************************************************************************/
/*
  Purpose:

    MONO_VALUE evaluates a monomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, unsigned M, the spatial dimension.

    Input, unsigned F[M], the exponents of the monomial.

    Input, double X[M*N], the coordinates of the evaluation points.

    Output, double MONO_VALUE[N], the value of the monomial at X.
*/
{
    unsigned m = _ID;
    unsigned i;
    unsigned char v = 1;

    for (i = 0; i < m; i++) {
        v = v * (unsigned char) pow(gf16v_get_ele(x, i), f[i]);//TODO: Casting not the best way?
        v = v % 16; //TODO: BULLSHIT, solve it with a multiply loop!
    }


    return v;
}

/******************************************************************************/

void perm_check0(unsigned n, const unsigned *p)

/******************************************************************************/
/*
  Purpose:

    PERM_CHECK0 checks a 0-based permutation.

  Discussion:

    The routine verifies that each of the ints from 0 to
    to N-1 occurs among the N entries of the permutation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, unsigned N, the number of entries.

    Input, unsigned P[N], the array to check.
*/
{
    unsigned ierror;
    unsigned location;
    unsigned value;

    for (value = 0; value < n; value++) {
        ierror = 1;

        for (location = 0; location < n; location++) {
            if (p[location] == value) {
                ierror = 0;
                break;
            }
        }

        if (ierror != 0) {
            fprintf(stderr, "\n");
            fprintf(stderr, "PERM_CHECK0 - Fatal error!\n");
            fprintf(stderr, "  Permutation is missing value %d\n", value);
            exit(1);
        }

    }
}

/******************************************************************************/

void
polynomial_add(unsigned char *destSummand, unsigned dest_offset, unsigned dest_grade, const unsigned char *summand,
               unsigned summand_offset, unsigned summand_o, const unsigned int *summand_e)

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_ADD adds two polynomials.

    - Attention: One summand-polynomial must not be the destination-polynomial at the same time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2013

  Author:

    John Burkardt

 4
  Parameters:

    Input, unsigned O1, the "order" of polynomial 1.

    Input, double C1[O1], the coefficients of polynomial 1.

    Input, unsigned E1[O1], the indices of the exponents of 
    polynomial 1.

    Input, unsigned O2, the "order" of polynomial 2.

    Input, double C2[O2], the coefficients of polynomial 2.

    Input, unsigned E2[O2], the indices of the exponents of 
    polynomial 2.

    Output, unsigned *O, the "order" of the polynomial sum.

    Output, double C[*O], the coefficients of the polynomial sum.

    Output, unsigned E[*O], the indices of the exponents of 
    the polynomial sum.
*/
{
    unsigned B_o = _grade_n_poly_terms(dest_grade);
    unsigned dest_o = summand_o + B_o;;

    unsigned dest_e[dest_o];

    //normal Polynom has standard e
    unsigned const *B_e = _full_e_power2;

    //needed, because else this function will write to far into the key..
    unsigned char tmp_B[(B_o + 1) / 2];
    unsigned char tmp_summand[(summand_o + 1) / 2];
    unsigned char tmp_dst[(summand_o + B_o + 1) / 2];

    gf16_copy(tmp_summand, 0, summand, summand_offset, summand_o);
    gf16_grade_n_poly_copy(tmp_B, 0, destSummand, dest_offset, dest_grade);

    polynomial_print(summand_o, tmp_summand, 0, summand_e, "summand:");
    polynomial_print(B_o, tmp_B, 0, _full_e_power2, "tmp_C:");

    /// MAIN PART///
    r8vec_concatenate(summand_o, tmp_summand, B_o, tmp_B, tmp_dst, 0);
    i4vec_concatenate(summand_o, summand_e, B_o, B_e, dest_e);

    polynomial_sort(dest_o, tmp_dst, 0, dest_e);

    polynomial_print(dest_o, tmp_dst, 0, dest_e, "sorted C:");

    polynomial_compress(dest_o, tmp_dst, 0, dest_e, &dest_o, tmp_dst, 0, dest_e);

    polynomial_print(dest_o, tmp_dst, 0, dest_e, "compressed C:");

    /// END ///

    gf16_copy(destSummand, dest_offset, tmp_dst, 0, dest_o);
}

/******************************************************************************/

void polynomial_compress(unsigned o1, unsigned char c1[], unsigned c1_offset, unsigned const e1[], unsigned *o2,
                         unsigned char c2[],
                         unsigned c2_offset,
                         unsigned e2[])

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_COMPRESS compresses a polynomial.

  Discussion:

    The function polynomial_sort ( ) should be called first, or else
    the E1 vector should be in ascending sorted order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, unsigned O1, the "order" of the polynomial.

    Input, double C1[O1], the coefficients of the polynomial.

    Input, unsigned E1[O1], the indices of the exponents of 
    the polynomial.

    Output, unsigned *O2, the "order" of the polynomial.

    Output, double C2[*O2], the coefficients of the polynomial.

    Output, unsigned E2[*O2], the indices of the exponents of 
    the polynomial.
*/
{
    unsigned char tmpA;
    unsigned char tmpB;
    unsigned char tmpSum;
    unsigned get;
    unsigned put;
/*
  Add coefficients associated with the same exponent.
*/
    get = 0;
    put = 0;

    while (get < o1) {
        get = get + 1;

        if (0 == put) {
            put = put + 1;
            gf16v_set_ele(c2, put - 1 + c2_offset, gf16v_get_ele(c1, get - 1 + c1_offset));
            e2[put - 1] = e1[get - 1];
        } else {
            if (e2[put - 1] == e1[get - 1]) {
                tmpA = gf16v_get_ele(c2, put - 1 + c2_offset);
                tmpB = gf16v_get_ele(c1, get - 1 + c1_offset);
                tmpSum = tmpA ^ tmpB;
                //TODO: check!!
                gf16v_set_ele(c2, put - 1, tmpSum);
            } else {
                put = put + 1;
                gf16v_set_ele(c2, put - 1 + c2_offset, gf16v_get_ele(c1, get - 1 + c1_offset));
                e2[put - 1] = e1[get - 1];
            }
        }
    }

    *o2 = put;
}

/******************************************************************************/

void
polynomial_mul(const unsigned char *factor_A, unsigned A_offset, unsigned A_grade, const unsigned char *factor_B,
               unsigned B_offset, unsigned B_grade, unsigned char *C, unsigned C_offset, unsigned *C_o,
               unsigned int C_e[])

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_MUL multiplies two polynomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, unsigned M, the spatial dimension. -> in this case same as _ID

    Input, unsigned O1, the "order" of polynomial 1.

    Input, double C1[O1], the coefficients of polynomial 1.

    Input, unsigned E1[O1], the indices of the exponents of 
    polynomial 1.

    Input, unsigned O2, the "order" of polynomial 2.

    Input, double C2[O2], the coefficients of polynomial 2.

    Input, unsigned E2[O2], the indices of the exponents of 
    polynomial 2.

    Output, unsigned *O, the "order" of the polynomial product.

    Output, double C[*O], the coefficients of the polynomial product.

    Output, unsigned E[*O], the indices of the exponents of the 
    polynomial product.
*/
{
    unsigned A_o = _grade_n_poly_terms(A_grade);
    unsigned const *A_e = _full_e_power2;

    unsigned B_o = _grade_n_poly_terms(B_grade);
    unsigned const *B_e = _full_e_power2;

    if (A_grade == 0) A_e = _lin_e_power2;
    if (B_grade == 0) B_e = _lin_e_power2;

    int m = _ID;

    unsigned f[m];
    unsigned f1[m];
    unsigned f2[m];
    unsigned i;
    unsigned j;
    int k;

    *C_o = 0;
    for (j = 0; j < B_o; j++) {
        for (i = 0; i < A_o; i++) {
            //c[*o] = c1[i] * c2[j]; -> adapt for GF16:

            unsigned char tmp_a = gf16v_get_ele(factor_A, i + A_offset);
            unsigned char tmp_b = gf16v_get_ele(factor_B, j + B_offset);
            unsigned char tmp_product = gf16_mul(tmp_a, tmp_b);

            gf16v_set_ele(C, *C_o + C_offset, tmp_product);

            mono_unrank_grlex(m, A_e[i], f1);
            mono_unrank_grlex(m, B_e[j], f2);
            for (k = 0; k < m; k++) {
                f[k] = f1[k] + f2[k];
            }
            C_e[*C_o] = mono_rank_grlex(m, f);
            memset(f1, 0, m * sizeof(unsigned));
            memset(f2, 0, m * sizeof(unsigned));
            *C_o = *C_o + 1;
        }
    }

    polynomial_sort(*C_o, C, C_offset, C_e);
    polynomial_compress(*C_o, C, C_offset, C_e, C_o, C, C_offset, C_e);
}

/******************************************************************************/

void polynomial_print(unsigned o, const unsigned char *c, unsigned gf16_offset, const unsigned *e, char *title)

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_PRint prints a polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, unsigned M, the spatial dimension.

    Input, unsigned O, the "order" of the polynomial, that is,
    simply the number of terms.

    Input, double C[O], the coefficients.

    Input, unsigned E[O], the indices of the exponents.

    Input, char *TITLE, a title.
*/
{
    int m = _ID;
    unsigned f[m];
    int i;
    unsigned j;

    printf("%s\n", title);


    if (o == 0) {
        printf("      0.\n");
    } else {
        for (j = 0; j < o; j++) {
            printf("    ");
            printf("+ ");
            printf("%hhu * x^(", gf16v_get_ele(c, j + gf16_offset));

            mono_unrank_grlex(m, e[j], f);
            for (i = 0; i < m; i++) {
                printf("%d", f[i]);
                if (i < m - 1) {
                    printf(",");
                } else {
                    printf(")");
                }
            }

            if (j == o - 1) {
                printf(".");
            }
            printf("\n");
        }
    }
}

/******************************************************************************/


void polynomial_sort(unsigned o, unsigned char c[], unsigned offset, unsigned e[])

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_SORT sorts the information in a polynomial.

  Discussion

    The coefficients C and exponents E are rearranged so that 
    the elements of E are in ascending order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, unsigned O, the "order" of the polynomial.

    Input/output, double C[O], the coefficients of the polynomial.

    Input/output, unsigned E[O], the indices of the exponents of 
    the polynomial.
*/
{
    unsigned *indx;

    indx = i4vec_sort_heap_index_a(o, e);

    i4vec_permute(o, indx, e); // indx=3
    r8vec_permute(o, indx, c, offset);

    free(indx);
}

/******************************************************************************/

unsigned char polynomial_value(unsigned o, const unsigned char *c, unsigned offset, unsigned const e[],
                               unsigned char *x)

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_VALUE evaluates a polynomial.

  Discussion:

    The polynomial is evaluated term by term, and no attempt is made to
    use an approach such as Horner's method to speed up the process.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, unsigned M, the spatial dimension.

    Input, unsigned O, the "order" of the polynomial.

    Input, double C[O], the coefficients of the polynomial.

    Input, unsigned E[O], the indices of the exponents 
    of the polynomial.

    Input, unsigned NX, the number of evaluation points.

    Input, double X[M*NX], the coordinates of the evaluation points.

    Output, double POLYNOMIAL_VALUE[NX], the value of the polynomial at X.
*/
{
    int m = _ID;
    unsigned f[m];
    unsigned j;
    unsigned char p;
    unsigned char v;

    p = 0;


    for (j = 0; j < o; j++) {
        mono_unrank_grlex(m, e[j], f);
        v = mono_value(f, x);
        p = p ^ gf16_mul(gf16v_get_ele(c, offset + j), v);
    }

    return p;
}

/******************************************************************************/

void
r8vec_concatenate(unsigned n1, const unsigned char a[], unsigned n2, const unsigned char b[], unsigned char c[],
                  unsigned offset)

/******************************************************************************/
/*
  Purpose:

    R8VEC_CONCATENATE concatenates two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, unsigned N1, the number of entries in the first vector.

    Input, double A[N1], the first vector.

    Input, unsigned N2, the number of entries in the second vector.

    Input, double B[N2], the second vector.

    Output, double C[N1+N2], the concatenated vector.
*/
{
    unsigned i;

    for (i = 0; i < n1; i++) {
        gf16v_set_ele(c, i + offset, gf16v_get_ele(a, i));
//        c[i] = a[i];
    }
    for (i = 0; i < n2; i++) {
        gf16v_set_ele(c, n1 + i + offset, gf16v_get_ele(b, i));
//        c[n1 + i] = b[i];
    }
}

/******************************************************************************/

void r8vec_permute(unsigned n1, unsigned p1[], unsigned char a[], unsigned offset)

/******************************************************************************/
/*
  Purpose:

    R8VEC_PERMUTE permutes an R8VEC in place.

  Discussion:

    An R8VEC is a vector of R8's.

    This routine permutes an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.

  Example:

    Input:

      N = 5
      P = (   1,   3,   4,   0,   2 )
      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )

    Output:

      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, unsigned N, the number of objects.

    Input, unsigned P[N], the permutation.

    Input/output, double A[N], the array to be permuted.
*/
{
    int *p = (int *) p1;
    int n = (int) n1;

    unsigned char a_temp;
    int i;
    int iget;
    int iput;
    int istart;

    perm_check0(n, p1);
/*
  In order for the sign negation trick to work, we need to assume that the
  entries of P are strictly positive.  Presumably, the lowest number is 0.
  So temporarily add 1 to each entry to force positivity.
*/
    for (i = 0; i < n; i++) {
        p[i] = p[i] + 1;
    }
/*
  Search for the next element of the permutation that has not been used.
*/
    for (istart = 1; istart <= n; istart++) {
        if (p[istart - 1] < 0) {
            continue;
        } else if (p[istart - 1] == istart) {
            p[istart - 1] = -p[istart - 1];
            continue;
        } else {
            a_temp = gf16v_get_ele(a, istart - 1 + offset);
            iget = istart;
/*
  Copy the new value into the vacated entry.
*/
            for (;;) {
                iput = iget;
                iget = p[iget - 1];

                p[iput - 1] = -p[iput - 1];

                if (iget < 1 || n < iget) {
                    fprintf(stderr, "\n");
                    fprintf(stderr, "R8VEC_PERMUTE - Fatal error!\n");
                    fprintf(stderr, "  A permutation index is out of range.\n");
                    fprintf(stderr, "  P(%d) = %d\n", iput, iget);
                    exit(1);
                }

                if (iget == istart) {
                    gf16v_set_ele(a, iput - 1 + offset, a_temp);
                    break;
                }
                gf16v_set_ele(a, iput - 1 + offset, gf16v_get_ele(a, iget - 1 + offset));
            }
        }
    }
/*
  Restore the signs of the entries.
*/
    for (i = 0; i < n; i++) {
        p[i] = -p[i];
    }
/*
  Restore the entries.
*/
    for (i = 0; i < n; i++) {
        p[i] = p[i] - 1;
    }
}