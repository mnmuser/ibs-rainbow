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


const unsigned _full_e_power2[N_QUARTIC_POLY(_ID)] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
const unsigned _lin_e_power2[_ID] = {2, 3};

/////////////////////////////////////////////////////////////////



static
void generate_S_T(unsigned char *s_and_t, prng_t *prng0) {
    prng_gen(prng0, s_and_t, _O1_BYTE * _O2 * _ID); // S1  => 16*32 = 512
    s_and_t += _O1_BYTE * _O2 * _ID;
    prng_gen(prng0, s_and_t, _V1_BYTE * _O1 * _ID); // T1
    s_and_t += _V1_BYTE * _O1 * _ID;
    prng_gen(prng0, s_and_t, _V1_BYTE * _O2 * _ID); // T2
    s_and_t += _V1_BYTE * _O2 * _ID;
    prng_gen(prng0, s_and_t, _O1_BYTE * _O2 * _ID); // T3
}


unsigned generate_l1_F12( unsigned char * sk, prng_t * prng0 ) {
    unsigned n_byte_generated = 0;
    prng_gen(prng0, sk, _O1_BYTE * N_TRIANGLE_TERMS(_V1) * _ID);
    /// Triangle: (n_var) (n_var*(n_var+1)/2)
    sk += _O1_BYTE * N_TRIANGLE_TERMS(_V1) * _ID;
    n_byte_generated += _O1_BYTE * N_TRIANGLE_TERMS(_V1) * _ID;

    prng_gen(prng0, sk, _O1_BYTE * _V1 * _O1 * _ID);  // l1_F2
    sk += _O1_BYTE * _V1 * _O1 * _ID;
    n_byte_generated += _O1_BYTE * _V1 * _O1 * _ID;

    return n_byte_generated;
}


static
unsigned generate_l2_F12356( unsigned char * sk, prng_t * prng0 ) {
    unsigned n_byte_generated = 0;

    prng_gen(prng0, sk, _O2_BYTE * N_TRIANGLE_TERMS(_V1) * _ID); // l2_F1
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_V1) * _ID;
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_V1) * _ID;

    prng_gen(prng0, sk, _O2_BYTE * _V1 * _O1 * _ID); // l2_F2
    sk += _O2_BYTE * _V1 * _O1 * _ID;
    n_byte_generated += _O2_BYTE * _V1 * _O1 * _ID;

    prng_gen(prng0, sk, _O2_BYTE * _V1 * _O2 * _ID); // l2_F3
    sk += _O2_BYTE * _V1 * _O1 * _ID;
    n_byte_generated += _O2_BYTE * _V1 * _O1 * _ID;

    prng_gen(prng0, sk, _O2_BYTE * N_TRIANGLE_TERMS(_O1) * _ID); // l2_F5
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_O1) * _ID;
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_O1) * _ID;

    prng_gen(prng0, sk, _O2_BYTE * _O1 * _O2 * _ID); // l2_F6
    n_byte_generated += _O2_BYTE * _O1 * _O2 * _ID;

    return n_byte_generated;
}

//TODO: Hier entsteht F
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
    for (unsigned i = 0; i < _O2; i++) {  /// t3 width
        gfmat_prod(temp, t1, _V1_BYTE, _O1, t3);
        gf256v_add(t4, temp, _V1_BYTE);
        t4 += _V1_BYTE;
        t3 += _O1_BYTE;
    }
}

static
void quartic_calculate_t4(unsigned char *t2_to_t4, const unsigned char *t1, const unsigned char *t3) {
    //  t4 = T_sk.t1 * T_sk.t3 - T_sk.t2
    unsigned char temp[(_O1_BYTE + 32) * N_QUADRATIC_POLY(_ID)];
    unsigned char *t4 = t2_to_t4;
    for (unsigned i = 0; i < _O2; i++) {  /// t3 width
        quartic_linear_gf16mat_prod_ref(temp, t1, _V1_BYTE, _O1, t3);
        //gf256v_add( t4 , temp , _V1_BYTE );

        unsigned e[25]; //e has to be long enough (?) for poly_add
        unsigned o = 0;
        unsigned char *tmp_t4 = malloc((N_QUADRATIC_POLY(_ID) + 1) / 2);
        for (unsigned j = 0; i < _O1; i++) {
            gf16_quadratic_poly_copy(tmp_t4, t4, j * _ID);
            polynomial_add(5, tmp_t4, _full_e_power2, 5, temp, _full_e_power2, &o, t4, 0, e);
        }
        free(tmp_t4);

        t4 += _V1_BYTE * N_QUADRATIC_POLY(_ID);
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
void quartic_obsfucate_l1_polys(unsigned char *l1_polys, const unsigned char *l2_polys, unsigned n_terms,
                                const unsigned char *s1) {
    unsigned char temp[_O1_BYTE * N_QUARTIC_POLY(_ID) + 32];
    while (n_terms--) { //for-loop *for* runaways
        quartic_gf16mat_prod_ref(temp, s1, _O1_BYTE, _O2, l2_polys); //(s1*l2 -> temp has grade4)
        unsigned e[25]; //e has to be long enough (?) for poly_add
        unsigned o = 0;
        unsigned char *tmp_l1_polys = malloc(N_QUARTIC_POLY((_ID) + 1) / 2);
        for (unsigned i = 0; i < _O1; i++) {
            gf16_quartic_poly_copy(tmp_l1_polys, 0, l1_polys, i * N_QUARTIC_POLY(_ID));
            polynomial_add(10, tmp_l1_polys, _full_e_power2, 15, temp, _full_e_power2, &o, l1_polys, 0, e);
        }
//        polynomial_print(o, l1_polys, 0, e, "temp:");
        free(tmp_l1_polys);
        //gf256v_add( l1_polys , temp , _O1_BYTE ); //add u32 (temp) on whole length of l1_polys
        l1_polys += _O1_BYTE * N_QUARTIC_POLY(_ID);
        l2_polys += _O2_BYTE * N_QUARTIC_POLY(_ID);
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


void generate_secretkey(msk_t *sk, const unsigned char *sk_seed) {
    _generate_secretkey(sk, sk_seed);
    calculate_t4(sk->t4, sk->t1, sk->t3);
}


void generate_keypair(mpk_t *rpk, msk_t *sk, const unsigned char *sk_seed) {

    _generate_secretkey(sk, sk_seed);

    // set up a temporary structure ext_cpk_t for calculating public key.
    ext_cpk_t *pk = (ext_cpk_t *) aligned_alloc(32, sizeof(ext_cpk_t));
    calculate_Q_from_F(pk, sk);

    /// at this point P = F o T ; S is still missing and P/Q is cubic on ID

    quartic_calculate_t4(sk->t4, sk->t1, sk->t3); // t4 = t4 + t1*t3

    obsfucate_l1_polys(pk->l1_Q1, pk->l2_Q1, N_TRIANGLE_TERMS(_V1), sk->s1); // -> integrate S :)
    quartic_obsfucate_l1_polys(pk->l1_Q2, pk->l2_Q2, _V1 * _O1, sk->s1);
    quartic_obsfucate_l1_polys(pk->l1_Q3, pk->l2_Q3, _V1 * _O2, sk->s1);
    quartic_obsfucate_l1_polys(pk->l1_Q5, pk->l2_Q5, N_TRIANGLE_TERMS(_O1), sk->s1);
    quartic_obsfucate_l1_polys(pk->l1_Q6, pk->l2_Q6, _O1 * _O2, sk->s1);
    quartic_obsfucate_l1_polys(pk->l1_Q9, pk->l2_Q9, N_TRIANGLE_TERMS(_O2), sk->s1);
    // so far, the pk contains the full pk but in ext_cpk_t format.

    quartic_extcpk_to_pk(rpk, pk);     // convert the public key from ext_cpk_t to mpk_t.

    free(pk);
}


////////////////////// IDENTITY ///////////////////////////

int calculate_usk(usk_t *usk, msk_t *msk, unsigned *id) {

}

int calculate_upk(upk_t *upk, mpk_t *mpk, unsigned char *id) {
    //TODO:
    calculate_values_public_key(upk,mpk,id);

    return 0;
}



void generate_identity_hash(unsigned char *digest, const unsigned char *id) {
    unsigned long long int id_length;
    id_length = sizeof(*id);
    hash_msg(digest, sizeof(*digest), id, id_length); // for simplicity I use the hash-function for messages
}

/// works because F,S,T and P are homogeneous
void multiply_identity_GF16(uint8_t *usk, uint8_t *upk, const unsigned char *id_hash, const uint8_t *msk,
                            const uint8_t *mpk) {
    int sk_size = sizeof(msk_t) - offsetof(msk_t, l1_F1);

    //Loop over sk
    multiply_ID_over_key(usk, msk, sk_size, id_hash);

    //Loop over pk
    multiply_ID_over_key(upk, mpk, CRYPTO_PUBLICKEYBYTES, id_hash);
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