/// @file rainbow_keypair_computation.h
/// @brief Functions for calculating pk/sk while generating keys.
///
/// Defining an internal structure of public key.
/// Functions for calculating pk/sk for key generation.
///

#ifndef _RAINBOW_KEYPAIR_COMP_H_
#define _RAINBOW_KEYPAIR_COMP_H_


#include "rainbow_keypair.h"


#ifdef  __cplusplus
extern  "C" {
#endif


/// @brief The (internal use) public key for rainbow
/// The (internal use) public key for rainbow. The public
/// polynomials are divided into l1_Q1, l1_Q2, ... l1_Q9,
/// l2_Q1, .... , l2_Q9.
///
typedef
struct rainbow_extend_publickey {
    unsigned char l1_Q1[_O1_BYTE * N_TRIANGLE_TERMS(_V1) * N_QUARTIC_POLY]; //TODO: necessary for normal PK!?
    unsigned char l1_Q2[_O1_BYTE * _V1 * _O1 * N_QUARTIC_POLY];
    unsigned char l1_Q3[_O1_BYTE * _V1 * _O2 * N_QUARTIC_POLY];
    unsigned char l1_Q5[_O1_BYTE * N_TRIANGLE_TERMS(_O1) * N_QUARTIC_POLY];
    unsigned char l1_Q6[_O1_BYTE * _O1 * _O2 * N_QUARTIC_POLY];
    unsigned char l1_Q9[_O1_BYTE * N_TRIANGLE_TERMS(_O2) * N_QUARTIC_POLY];
    /// -> 372480
    unsigned char l2_Q1[_O2_BYTE * N_TRIANGLE_TERMS(_V1) * N_QUARTIC_POLY];
    unsigned char l2_Q2[_O2_BYTE * _V1 * _O1 * N_QUARTIC_POLY];
    unsigned char l2_Q3[_O2_BYTE * _V1 * _O2 * N_QUARTIC_POLY];
    unsigned char l2_Q5[_O2_BYTE * N_TRIANGLE_TERMS(_O1) * N_QUARTIC_POLY];
    unsigned char l2_Q6[_O2_BYTE * _O1 * _O2 * N_QUARTIC_POLY];
    unsigned char l2_Q9[_O2_BYTE * N_TRIANGLE_TERMS(_O2) * N_QUARTIC_POLY];
} ext_cpk_t; /// TODO: normal: 148992, same as PK



///
/// @brief converting formats of public keys : from ext_cpk_t version to mpk_t
///
/// @param[out] pk       - the classic public key.
/// @param[in]  cpk      - the internal public key.
///
void extcpk_to_pk(mpk_t *pk, const ext_cpk_t *cpk);

void quartic_extcpk_to_pk(mpk_t *pk, const ext_cpk_t *cpk);

/////////////////////////////////////////////////

///
/// @brief Computing public key from secret key, Q = F o T
///
/// @param[out] Qs       - the public key
/// @param[in]  Fs       - the secret key
///
void calculate_Q_from_F(ext_cpk_t *Qs, const msk_t *Fs);


///
/// @brief Computing parts of the sk from parts of pk and sk
///
/// @param[out] Fs       - parts of the sk: l1_F1, l1_F2, l2_F1, l2_F2, l2_F3, l2_F5, l2_F6
/// @param[in]  Qs       - parts of the pk: l1_Q1, l1_Q2, l2_Q1, l2_Q2, l2_Q3, l2_Q5, l2_Q6
/// @param[in]  Ts       - parts of the sk: T1, T4, T3
///
void calculate_F_from_Q(msk_t *Fs, const msk_t *Qs, msk_t *Ts);

void write_gf16_to_quartic(unsigned char *q, const unsigned char *f, unsigned long length_f);

void set_quartic_zero(unsigned char *q, unsigned length);

void gf16_grade_n_poly_copy(unsigned char *dest, unsigned gf16_offset_dest, const unsigned char *src,
                            unsigned gf16_offset_src, unsigned grade_n);

#ifdef  __cplusplus
}
#endif

#endif  // _RAINBOW_KEYPAIR_COMP_H_

