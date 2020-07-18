/// @file rainbow_keypair.h
/// @brief Formats of key pairs and functions for generating key pairs.
/// Formats of key pairs and functions for generating key pairs.
///

#ifndef _RAINBOW_KEYPAIR_H_
#define _RAINBOW_KEYPAIR_H_


#include "rainbow_config.h"
#include <stdint.h>

/// Number of quadratic polynomials in n_var / Gaußsche Summenformel, auch zur Berechnung von Dreieckszahlen
#define N_TRIANGLE_TERMS(n_var) (n_var*(n_var+1)/2)

/// number of coefficients needed in public key (Quartic polynomial) -> Ding p.1143
#define N_QUARTIC_POLY ((_ID+1)*(_ID+2)*(_ID+3)*(_ID+4)/24)
/// number of terms in quadratic polynom
#define N_QUADRATIC_POLY ((_ID+1)*(_ID+2)/2)

#define N_CUBIC_POLY ((_ID+1)*(_ID+2)*(_ID+3)/6)

#define N_LINEAR_POLY ((_ID+1))

#ifdef  __cplusplus
extern  "C" {
#endif

/// standard e TODO: only for ID2


extern const unsigned _full_e_power2[N_QUARTIC_POLY + 3];

extern const unsigned _lin_e_power2[_ID];

unsigned _grade_n_poly_terms(unsigned grade); //TODO: extern, static, global namespace...

/// @brief master public key for classic rainbow
///
///  master public key for classic rainbow
///
typedef
struct rainbow_master_publickey {
    unsigned char pk[(_PUB_M_BYTE) * N_TRIANGLE_TERMS(_PUB_N) * N_QUARTIC_POLY]; // -> seems legit
    /// _Pub_M_Byte = (O1+O2)/2 = 32 -> Anzahl der Gleichungen (O1+O2) durch 2 wg GF16
    /// _ PUB_N = V1+O1+O2 = 32+32+32= 96 ->
    /// ((O1+O2)/2) * ((V1+O1+O2)*(V1+O1+O2+1)/2)
    /// 32 * (96 * 97 / 2)
    /// =unsigned char[148992] [CHECKED]
    /// * ID_len in QUAT_POLY
    /// p.16 in documentation
    /// same length as extend_publickey
} mpk_t;


/// @brief master secret key for classic rainbow
///
/// master secret key for classic rainbow
///
typedef
struct rainbow_master_secretkey {
    /// * ID_LEN...
    /// seed for generating secret key.
    /// Generating S, T, and F for classic rainbow.
    /// Generating S and T only for cyclic rainbow.
    unsigned char sk_seed[LEN_SKSEED];

    unsigned char s1[_O1_BYTE * _O2 * _ID];   ///< part of S map (upper right corner in matrix, p.3 in doc)
    unsigned char t1[_V1_BYTE * _O1 * _ID];   ///< part of T map
    unsigned char t4[_V1_BYTE * _O2 * _ID];   ///< part of T map (gets quadratic) -> also known as t2
    unsigned char t3[_O1_BYTE * _O2 * _ID];   ///< part of T map

    unsigned char l1_F1[_O1_BYTE * N_TRIANGLE_TERMS(_V1) * _ID];  ///< part of C-map, F1, Layer1
    unsigned char l1_F2[_O1_BYTE * _V1 * _O1 * _ID];                ///< part of C-map, F2, Layer1

    unsigned char l2_F1[_O2_BYTE * N_TRIANGLE_TERMS(_V1) * _ID];  ///< part of C-map, F1, Layer2
    unsigned char l2_F2[_O2_BYTE * _V1 * _O1 * _ID];                ///< part of C-map, F2, Layer2

    unsigned char l2_F3[_O2_BYTE * _V1 * _O2 * _ID];                ///< part of C-map, F3, Layer2
    unsigned char l2_F5[_O2_BYTE * N_TRIANGLE_TERMS(_O1) * _ID];  ///< part of C-map, F5, Layer2
    unsigned char l2_F6[_O2_BYTE * _O1 * _O2 * _ID];                ///< part of C-map, F6, Layer2
} msk_t;


/// @brief user public key for classic rainbow
///
///  user public key for classic rainbow
///
typedef
struct rainbow_user_publickey {
    unsigned char pk[(_PUB_M_BYTE) * N_TRIANGLE_TERMS(_PUB_N)];
} upk_t;


/// @brief user secret key for classic rainbow
///
/// user secret key for classic rainbow
///
typedef
struct rainbow_user_secretkey {
    ///
    /// seed for generating secret key.
    /// Generating S, T, and F for classic rainbow.
    /// Generating S and T only for cyclic rainbow.
    unsigned char sk_seed[LEN_SKSEED]; //TODO: hust?

    unsigned char s1[_O1_BYTE * _O2];   ///< part of S map
    unsigned char t1[_V1_BYTE * _O1];   ///< part of T map
    unsigned char t4[_V1_BYTE * _O2];   ///< part of T map
    unsigned char t3[_O1_BYTE * _O2];   ///< part of T map

    unsigned char l1_F1[_O1_BYTE * N_TRIANGLE_TERMS(_V1)];  ///< part of C-map, F1, Layer1
    unsigned char l1_F2[_O1_BYTE * _V1 * _O1];                ///< part of C-map, F2, Layer1

    unsigned char l2_F1[_O2_BYTE * N_TRIANGLE_TERMS(_V1)];  ///< part of C-map, F1, Layer2
    unsigned char l2_F2[_O2_BYTE * _V1 * _O1];                ///< part of C-map, F2, Layer2

    unsigned char l2_F3[_O2_BYTE * _V1 * _O2];                ///< part of C-map, F3, Layer2
    unsigned char l2_F5[_O2_BYTE * N_TRIANGLE_TERMS(_O1)];  ///< part of C-map, F5, Layer2
    unsigned char l2_F6[_O2_BYTE * _O1 * _O2];                ///< part of C-map, F6, Layer2
} usk_t;



/////////////////////////////////////


///
/// @brief Generate key pairs for classic rainbow.
///
/// @param[out] pk        - the public key.
/// @param[out] sk        - the secret key.
/// @param[in]  sk_seed   - seed for generating the secret key.
///
void generate_keypair(mpk_t *pk, msk_t *sk, const unsigned char *sk_seed);


////////////////////////////////////

///
/// @brief Generate secret key for classic rainbow.
///
/// @param[out] sk        - the secret key.
/// @param[in]  sk_seed   - seed for generating the secret key.
///
void generate_secretkey(msk_t *sk, const unsigned char *sk_seed);


////////////////////////////////////

int calculate_usk(usk_t *usk, msk_t *msk, unsigned char *id);

int calculate_upk(upk_t *upk, mpk_t *mpk, unsigned char *id);

/////////////////ID///////////////////

void generate_identity_hash(unsigned char *digest, const unsigned char *id, unsigned id_length);

void multiply_identity_sk(msk_t *usk, const unsigned char *id_hash, msk_t *msk);

/// Function for multiplying ID bitwise over the whole keys to get USK and UPK in GF16
/// \param usk
/// \param upk
/// \param id_hash
/// \param msk
/// \param mpk
void multiply_identity_GF16(uint8_t *usk, uint8_t *upk, const unsigned char *id_hash, const uint8_t *msk,
                            const uint8_t *mpk);


void
multiply_ID_over_key(uint8_t *dest_key, const uint8_t *key, unsigned long key_length, const unsigned char *id_hash);

#ifdef  __cplusplus
}
#endif

#endif //  _RAINBOW_KEYPAIR_H_
