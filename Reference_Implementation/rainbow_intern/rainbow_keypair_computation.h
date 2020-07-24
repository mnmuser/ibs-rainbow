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
struct rainbow_extend_masterpublickey {
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
} ext_mpk_t; /// 148992 * Quartic, same as PK

typedef
struct rainbow_extend_publickey {
    unsigned char l1_Q1[_O1_BYTE * N_TRIANGLE_TERMS(_V1)];
    unsigned char l1_Q2[_O1_BYTE * _V1 * _O1];
    unsigned char l1_Q3[_O1_BYTE * _V1 * _O2];
    unsigned char l1_Q5[_O1_BYTE * N_TRIANGLE_TERMS(_O1)];
    unsigned char l1_Q6[_O1_BYTE * _O1 * _O2];
    unsigned char l1_Q9[_O1_BYTE * N_TRIANGLE_TERMS(_O2)];

    unsigned char l2_Q1[_O2_BYTE * N_TRIANGLE_TERMS(_V1)];
    unsigned char l2_Q2[_O2_BYTE * _V1 * _O1];
    unsigned char l2_Q3[_O2_BYTE * _V1 * _O2];
    unsigned char l2_Q5[_O2_BYTE * N_TRIANGLE_TERMS(_O1)];
    unsigned char l2_Q6[_O2_BYTE * _O1 * _O2];
    unsigned char l2_Q9[_O2_BYTE * N_TRIANGLE_TERMS(_O2)];
} ext_cpk_t; /// size: unsigned char[148992]



///
/// @brief converting formats of public keys : from ext_mpk_t version to mpk_t
///
/// @param[out] pk       - the classic public key.
/// @param[in]  cpk      - the internal public key.
///
void extcpk_to_pk(upk_t *pk, const ext_cpk_t *cpk);

void quartic_extcpk_to_pk(mpk_t *pk, const ext_mpk_t *cpk);

/////////////////////////////////////////////////

///
/// @brief Computing public key from secret key, Q = F o T
///
/// @param[out] Qs       - the public key
/// @param[in]  Fs       - the secret key
///
void quartic_calculate_Q_from_F(ext_mpk_t *Qs, const msk_t *Fs);

void calculate_Q_from_F(ext_cpk_t *Qs, const usk_t *Fs, const usk_t *Ts);

/// @brief function to write from a linear funtion to a quartic function respecting gaps
/// \param q
/// \param f
/// \param length_f
void write_lin_to_quartic(unsigned char *q, const unsigned char *f, unsigned long length_f);

/// @brief function to write from a linear function without constant to a quadratic function respecting gaps
/// \param dest
/// \param src
/// \param length_f
void write_lin_wo_const_to_quadratic(unsigned char *dest, const unsigned char *src, unsigned long length_f);

void set_quartic_zero(unsigned char *q, unsigned length);

void gf16_grade_n_poly_copy(unsigned char *dest, unsigned gf16_offset_dest, const unsigned char *src,
                            unsigned gf16_offset_src, unsigned grade_n);

void gf16_copy(unsigned char *dest, unsigned gf16_offset_dest, const unsigned char *src,
               unsigned gf16_offset_src, unsigned size_values);

#ifdef  __cplusplus
}
#endif

#endif  // _RAINBOW_KEYPAIR_COMP_H_

