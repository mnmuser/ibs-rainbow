/// @file rainbow_keypair.c
/// @brief implementations of functions in rainbow_keypair.h
///
#include "rainbow_keypair.h"
#include "rainbow_keypair_computation.h"

#include "blas_comm.h"
#include "blas.h"
#include "rainbow_blas.h"

#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include "utils_prng.h"
#include "utils_hash.h"

#include "api.h"
#include "polynomial.h"


const unsigned _full_e[15] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
const unsigned _lin_e[2] = {2, 3};

unsigned _grade_n_poly_terms(unsigned grade) {
    switch (grade) {
        case 0:
            return _ID; //this is mathematically viewn not correct, but needed in case of linear polynom without constant factor (aka sk)
        case 1:
            return N_LINEAR_POLY;
        case 2:
            return N_QUADRATIC_POLY;
        case 3:
            return N_CUBIC_POLY;
        case 4:
            return N_QUARTIC_POLY;
        default:
            return 0;
    }
}


/////////////////////////////////////////////////////////////////



static
void generate_S_T(unsigned char *s_and_t, prng_t *prng0) {
    prng_gen(prng0, s_and_t, _O1_BYTE * _O2 * _ID); // S1  => 16*32 = 512
    s_and_t += _O1_BYTE * _O2 * _ID;
    prng_gen(prng0, s_and_t, _V1_BYTE * _O1 * _ID); // T1
    s_and_t += _V1_BYTE * _O1 * _ID;
    prng_gen(prng0, s_and_t, _V1_BYTE * _O2 * _ID); // T2 bzw. T4
    s_and_t += _V1_BYTE * _O2 * _ID; // skip more fields (needed for later (t4-calculation))
    prng_gen(prng0, s_and_t, _O1_BYTE * _O2 * _ID); // T3
}


unsigned generate_l1_F12( unsigned char * sk, prng_t * prng0 ) {
    unsigned n_byte_generated = 0;
    prng_gen(prng0, sk, _O1_BYTE * N_TRIANGLE_TERMS(_V1) * N_LINEAR_POLY);
    /// Triangle: (n_var) (n_var*(n_var+1)/2)
    sk += _O1_BYTE * N_TRIANGLE_TERMS(_V1) * N_LINEAR_POLY;
    n_byte_generated += _O1_BYTE * N_TRIANGLE_TERMS(_V1) * N_LINEAR_POLY;

    prng_gen(prng0, sk, _O1_BYTE * _V1 * _O1 * N_LINEAR_POLY);  // l1_F2
    sk += _O1_BYTE * _V1 * _O1 * N_LINEAR_POLY;
    n_byte_generated += _O1_BYTE * _V1 * _O1 * N_LINEAR_POLY;

    return n_byte_generated;
}


static
unsigned generate_l2_F12356( unsigned char * sk, prng_t * prng0 ) {
    unsigned n_byte_generated = 0;

    prng_gen(prng0, sk, _O2_BYTE * N_TRIANGLE_TERMS(_V1) * N_LINEAR_POLY); // l2_F1
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_V1) * N_LINEAR_POLY;
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_V1) * N_LINEAR_POLY;

    prng_gen(prng0, sk, _O2_BYTE * _V1 * _O1 * N_LINEAR_POLY); // l2_F2
    sk += _O2_BYTE * _V1 * _O1 * N_LINEAR_POLY;
    n_byte_generated += _O2_BYTE * _V1 * _O1 * N_LINEAR_POLY;

    prng_gen(prng0, sk, _O2_BYTE * _V1 * _O2 * N_LINEAR_POLY); // l2_F3
    sk += _O2_BYTE * _V1 * _O1 * N_LINEAR_POLY;
    n_byte_generated += _O2_BYTE * _V1 * _O1 * N_LINEAR_POLY;

    prng_gen(prng0, sk, _O2_BYTE * N_TRIANGLE_TERMS(_O1) * N_LINEAR_POLY); // l2_F5
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_O1) * N_LINEAR_POLY;
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_O1) * N_LINEAR_POLY;

    prng_gen(prng0, sk, _O2_BYTE * _O1 * _O2 * N_LINEAR_POLY); // l2_F6
    n_byte_generated += _O2_BYTE * _O1 * _O2 * N_LINEAR_POLY;

    return n_byte_generated;
}

// Hier entsteht F
static
void generate_B1_B2( unsigned char * sk , prng_t * prng0 ) {
    sk += generate_l1_F12(sk, prng0); //Layer 1
    generate_l2_F12356(sk, prng0); // Layer 2
}


/////////////////////////////////////////////////////////

static
void calculate_t4( unsigned char * t2_to_t4 , const unsigned char *t1 , const unsigned char *t3 )
{
    //  t4 = T_sk.t1 * T_sk.t3 - T_sk.t2
    unsigned char temp[_V1_BYTE + 32];
    unsigned char *t4 = t2_to_t4;
    for (unsigned i = 0; i < _O2; i++) {  // t3 width
        gfmat_prod(temp, t1, _V1_BYTE, _O1, t3);
        gf256v_add(t4, temp, _V1_BYTE);
        t4 += _V1_BYTE;
        t3 += _O1_BYTE;
    }
}

static
void quadratic_calculate_t4(unsigned char *t2_to_t4, const unsigned char *t1, const unsigned char *t3) {
    unsigned char tmp_t2[_V1_BYTE * _O2 * _ID];
    memcpy(tmp_t2, t2_to_t4, _V1_BYTE * _O2 * _ID);
    write_lin_wo_const_to_quadratic(t2_to_t4, tmp_t2, _V1_BYTE * _O2);

    //  t4 = T_sk.t1 * T_sk.t3 - T_sk.t2
    unsigned char temp[_V1_BYTE * N_QUADRATIC_POLY + 32];
    unsigned char *t4 = t2_to_t4;
    for (unsigned i = 0; i < _V1_BYTE; i++) {
        quadratic_gf16mat_prod_ref(temp, t1, _V1_BYTE, _O1, t3);
        //gfmat_prod(temp, t1, _V1_BYTE, _O1, t3);
        for (unsigned i = 0; i < _O1_BYTE * 2; i++) {
            polynomial_add(t2_to_t4, i * N_QUADRATIC_POLY, 1, temp, i * _ID, N_QUADRATIC_POLY, _full_e);
        }
//        gf256v_add(t4, temp, _V1_BYTE);
        t4 += _V1_BYTE * N_QUADRATIC_POLY;
        t3 += _O1_BYTE * _ID;
    }
}

static
void
obsfucate_l1_polys(unsigned char *l1_polys, const unsigned char *l2_polys, unsigned n_terms, const unsigned char *s1) {
    unsigned char temp[_O1_BYTE + 32];
    while (n_terms--) {
        gfmat_prod(temp, s1, _O1_BYTE, _O2, l2_polys);
        gf256v_add(l1_polys, temp, _O1_BYTE);
        l1_polys += _O1_BYTE;
        l2_polys += _O2_BYTE;
    }
}

static
void quartic_obsfucate_l1_polys(unsigned char *l1_polys, const unsigned char *l2_polys, unsigned poly_grade,
                                unsigned n_terms,
                                const unsigned char *s1) {
    unsigned char temp[_O1_BYTE * N_QUARTIC_POLY + 32];

    while (n_terms--) { //for-loop *for* runaways
        quartic_gf16mat_prod_ref(temp, s1, _O1_BYTE, _O2, l2_polys, poly_grade);

        for (unsigned i = 0; i < _O1_BYTE * 2; i++) {
            polynomial_add(l1_polys, i * N_QUARTIC_POLY, poly_grade, temp, i * N_QUARTIC_POLY, N_QUARTIC_POLY,
                           _full_e);
        }
        //gf256v_add( l1_polys , temp , _O1_BYTE ); //add u32 (temp) on whole length of l1_polys

        l1_polys += _O1_BYTE * N_QUARTIC_POLY;
        l2_polys += _O2_BYTE * N_QUARTIC_POLY;
    }
}


///////////////////  Classic //////////////////////////////////


static
void _generate_secretkey(msk_t *sk, const unsigned char *sk_seed) {
    memcpy(sk->sk_seed, sk_seed, LEN_SKSEED);

    // set up prng
    prng_t prng0;
    prng_set(&prng0, sk_seed, LEN_SKSEED);

    // generating secret key with prng.
    generate_S_T(sk->s1, &prng0);
    generate_B1_B2(sk->l1_F1, &prng0);

    // clean prng
    memset( &prng0 , 0 , sizeof(prng_t) );
}


void generate_keypair(mpk_t *rpk, msk_t *sk, const unsigned char *sk_seed) {

    _generate_secretkey(sk, sk_seed);

    // set up a temporary structure ext_mpk_t for calculating public key.
    ext_mpk_t *pk = (ext_mpk_t *) aligned_alloc(32, sizeof(ext_mpk_t));
    quartic_calculate_Q_from_F(pk, sk);

    /// at this point P = F o T ; S is still missing and P/Q is cubic on ID

    ///TODO: is this right? (let's test)
//    quadratic_calculate_t4(sk->t4, sk->t1, sk->t3); // t2 = t2 + t1*t3

    quartic_obsfucate_l1_polys(pk->l1_Q1, pk->l2_Q1, 1, N_TRIANGLE_TERMS(_V1), sk->s1); // -> integrate S :)

    ///CHECK Q1
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q1, 0, _full_e, "obsfucated L1_Q1 ");
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q1, sizeof(pk->l1_Q1) * 2 - N_QUARTIC_POLY, _full_e,
                     "obsfucated L1_Q1 end");
    ///

    quartic_obsfucate_l1_polys(pk->l1_Q2, pk->l2_Q2, 2, _V1 * _O1, sk->s1);

    ///CHECK Q2
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q2, 0, _full_e, "obsfucated L1_Q2: ");
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q2, sizeof(pk->l1_Q2) * 2 - N_QUARTIC_POLY, _full_e,
                     "obsfucated L1_Q2 end: "); //last position in Q2
    ///

    quartic_obsfucate_l1_polys(pk->l1_Q3, pk->l2_Q3, 2, _V1 * _O2, sk->s1);

    ///CHECK Q3
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q3, 0, _full_e, "obsfucated Q3:");
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q3, sizeof(pk->l1_Q3) * 2 - N_QUARTIC_POLY, _full_e,
                     "obsfucated L1_Q3 end: ");
    ///

    quartic_obsfucate_l1_polys(pk->l1_Q5, pk->l2_Q5, 3, N_TRIANGLE_TERMS(_O1), sk->s1);

    ///CHECK Q5
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q5, 0, _full_e, "obsfucated l1_Q5(0): ");
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q5, sizeof(pk->l1_Q5) * 2 - N_QUARTIC_POLY, _full_e,
                     "obsfucated L1_Q5 end");
    ///

    quartic_obsfucate_l1_polys(pk->l1_Q6, pk->l2_Q6, 3, _O1 * _O2, sk->s1);

    ///CHECK Q6
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q6, 0, _full_e, "obsfucated Q6:");
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q6, sizeof(pk->l1_Q6) * 2 - N_QUARTIC_POLY, _full_e,
                     "obsfucated L1_Q6 end: ");
    ///

    quartic_obsfucate_l1_polys(pk->l1_Q9, pk->l2_Q9, 3, N_TRIANGLE_TERMS(_O2), sk->s1);

    ///CHECK Q9
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q9, 0, _full_e, "obsfucated l1_Q9(0): ");
    polynomial_print(N_QUARTIC_POLY, pk->l1_Q9, sizeof(pk->l1_Q9) * 2 - N_QUARTIC_POLY, _full_e,
                     "obsfucated L1_Q9 end");
    ///

    // so far, the pk contains the full pk but in ext_mpk_t format.

    quartic_extcpk_to_pk(rpk, pk);   //TODO comment in  // convert the public key from ext_mpk_t to mpk_t.

    free(pk);
}


////////////////////// IDENTITY ///////////////////////////

int calculate_usk(usk_t *usk, msk_t *msk, unsigned char *id) {
    calculate_values_secret_key((unsigned char *) usk, (unsigned char *) msk, id);

    //last but not least: (so I don't need to make t4 quartic)
    calculate_t4(usk->t4, usk->t1, usk->t3);
    return 0;
}

int calculate_upk(upk_t *upk, mpk_t *mpk, unsigned char *id) {
    calculate_values_public_key((unsigned char *) upk, (unsigned char *) mpk, id);
    return 0;
}

/// for debugging purposes
/// \param upk
/// \param usk
/// \return
int calculate_usk_and_upk(usk_t *usk, upk_t *upk, msk_t *msk, unsigned char *id) {
    //calculate usk
    calculate_values_secret_key((unsigned char *) usk, (unsigned char *) msk, id);

    //calculate upk
    ext_cpk_t *pk = (ext_cpk_t *) aligned_alloc(32, sizeof(ext_cpk_t));
    calculate_Q_from_F(pk, usk, usk);
    obsfucate_l1_polys(pk->l1_Q1, pk->l2_Q1, N_TRIANGLE_TERMS(_V1), usk->s1);
    obsfucate_l1_polys(pk->l1_Q2, pk->l2_Q2, _V1 * _O1, usk->s1);
    obsfucate_l1_polys(pk->l1_Q3, pk->l2_Q3, _V1 * _O2, usk->s1);
    obsfucate_l1_polys(pk->l1_Q5, pk->l2_Q5, N_TRIANGLE_TERMS(_O1), usk->s1);
    obsfucate_l1_polys(pk->l1_Q6, pk->l2_Q6, _O1 * _O2, usk->s1);
    obsfucate_l1_polys(pk->l1_Q9, pk->l2_Q9, N_TRIANGLE_TERMS(_O2), usk->s1);
    extcpk_to_pk(upk, pk);
//    memcpy(upk,pk,sizeof(upk_t));
    free(pk);

    //t4 in usk
    calculate_t4(usk->t4, usk->t1, usk->t3);
    return 0;
}

void generate_identity_hash(unsigned char *digest, unsigned char *id, unsigned id_length) {
    unsigned hash_length = (_ID + 1) / 2;
    hash_msg(digest, hash_length, id, id_length); // for simplicity I use the hash-function for messages
}

/// works because F,S,T and P are homogeneous
void multiply_identity_GF16(uint8_t *usk, uint8_t *upk, const unsigned char *id_hash, const uint8_t *msk,
                            const uint8_t *mpk) {
    unsigned sk_size = sizeof(msk_t) - offsetof(msk_t, l1_F1);

    //Loop over sk
    multiply_ID_over_key(usk, msk, sk_size, id_hash);

    //Loop over pk
    multiply_ID_over_key(upk, mpk, CRYPTO_MASTER_PUBLIC_KEY_BYTES, id_hash);
}

void multiply_ID_over_key(unsigned char *dest_key, const unsigned char *key, const unsigned long key_length,
                          const unsigned char *id_hash) {
    unsigned long length_in_elements = key_length * 2;

    for (unsigned i = 0; i < length_in_elements; i++) {
        //get element from key at i
        uint8_t key_element = gf16v_get_ele(key, i);
        //multiply them
        /// With that way we only have 16 identities
        /// XOR is not working
        /// addition is not working
        /// corresponding element out of ID-Hash is not working
        /// multiply with a integer in GF16 works
        uint8_t product = gf16_mul(key_element, gf16v_get_ele(id_hash, 0));
        //set them for the user key
        gf16v_set_ele(dest_key, i, product);
    }
}

static const unsigned *fill_full_e() {
    static unsigned e[N_QUARTIC_POLY]; //static, can be used everywhere without the risk of malloc/free
    for (unsigned i = 0; i < N_QUARTIC_POLY; ++i) {
        e[i] = i + 1;
    }
    return e;
}