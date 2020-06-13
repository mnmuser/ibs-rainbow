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
#define N_QUARTIC_POLY(id_var) ((id_var+1)*(id_var+2)*(id_var+3)*(id_var+4)/24)
/// number of terms in quadratic polynom
#define N_QUADRATIC_POLY(id_var) ((id_var+1)*(id_var+2)/2)

#define N_CUBIC_POLY(id_var) ((id_var+1)*(id_var+2)*(id_var+3)/6)

#ifdef  __cplusplus
extern  "C" {
#endif

/// standard e TODO: only for ID2
const unsigned full_e_power2[N_QUARTIC_POLY(_ID)] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

/// @brief public key for classic rainbow
///
///  public key for classic rainbow
///
typedef
struct rainbow_publickey {
    unsigned char pk[(_PUB_M_BYTE) * N_TRIANGLE_TERMS(_PUB_N) * N_QUARTIC_POLY(_ID)]; // -> seems legit
    /// _Pub_M_Byte = (O1+O2)/2 = 32 -> Anzahl der Gleichungen (O1+O2) durch 2 wg GF16
    /// _ PUB_N = V1+O1+O2 = 32+32+32= 96 ->
    /// ((O1+O2)/2) * ((V1+O1+O2)*(V1+O1+O2+1)/2)
    /// 32 * (96 * 97 / 2)
    /// =unsigned char[148992] [CHECKED]
    /// * ID_len in QUAT_POLY
    /// p.16 in documentation
    /// same length as extend_publickey
} pk_t;


/// @brief secret key for classic rainbow
///
/// secret key for classic rainbow
///
typedef
struct rainbow_secretkey {
    /// * ID_LEN...
    /// seed for generating secret key.
    /// Generating S, T, and F for classic rainbow.
    /// Generating S and T only for cyclic rainbow.
    unsigned char sk_seed[LEN_SKSEED];

    unsigned char s1[_O1_BYTE * _O2 * _ID];   ///< part of S map (upper right corner in matrix, p.3 in doc)
    unsigned char t1[_V1_BYTE * _O1 * _ID];   ///< part of T map
    unsigned char t4[_V1_BYTE * _O2 * N_QUADRATIC_POLY(_ID)];   ///< part of T map (gets quadratic)
    unsigned char t3[_O1_BYTE * _O2 * _ID];   ///< part of T map

    unsigned char l1_F1[_O1_BYTE * N_TRIANGLE_TERMS(_V1) * _ID];  ///< part of C-map, F1, Layer1
    unsigned char l1_F2[_O1_BYTE * _V1 * _O1 * _ID];                ///< part of C-map, F2, Layer1

    unsigned char l2_F1[_O2_BYTE * N_TRIANGLE_TERMS(_V1) * _ID];  ///< part of C-map, F1, Layer2
    unsigned char l2_F2[_O2_BYTE * _V1 * _O1 * _ID];                ///< part of C-map, F2, Layer2

    unsigned char l2_F3[_O2_BYTE * _V1 * _O2 * _ID];                ///< part of C-map, F3, Layer2
    unsigned char l2_F5[_O2_BYTE * N_TRIANGLE_TERMS(_O1) * _ID];  ///< part of C-map, F5, Layer2
    unsigned char l2_F6[_O2_BYTE * _O1 * _O2 * _ID];                ///< part of C-map, F6, Layer2
} sk_t;


////////// CYCLIC ///////////////////////

/// @brief public key for cyclic rainbow
///
///  public key for cyclic rainbow
///
typedef
struct rainbow_publickey_cyclic {
    unsigned char pk_seed[LEN_PKSEED];                      ///< seed for generating l1_Q1,l1_Q2,l2_Q1,l2_Q2,l2_Q3,l2_Q5,l2_Q6

    unsigned char l1_Q3[_O1_BYTE * _V1 * _O2];                ///< Q3, layer1
    unsigned char l1_Q5[_O1_BYTE * N_TRIANGLE_TERMS(_O1)];  ///< Q5, layer1
    unsigned char l1_Q6[_O1_BYTE * _O1 * _O2];                ///< Q6, layer1
    unsigned char l1_Q9[_O1_BYTE * N_TRIANGLE_TERMS(_O2)];  ///< Q9, layer1

    unsigned char l2_Q9[_O2_BYTE * N_TRIANGLE_TERMS(_O2)];  ///< Q9, layer2
} cpk_t;



/// @brief compressed secret key for cyclic rainbow
///
/// compressed secret key for cyclic rainbow
///
typedef
struct rainbow_secretkey_cyclic {
    unsigned char pk_seed[LEN_PKSEED];   ///< seed for generating a part of public key.
    unsigned char sk_seed[LEN_SKSEED];   ///< seed for generating a part of secret key.
} csk_t;



/////////////////////////////////////


///
/// @brief Generate key pairs for classic rainbow.
///
/// @param[out] pk        - the public key.
/// @param[out] sk        - the secret key.
/// @param[in]  sk_seed   - seed for generating the secret key.
///
void generate_keypair(pk_t *pk, sk_t *sk, const unsigned char *sk_seed);

///
/// @brief Generate key pairs for cyclic rainbow.
///
/// @param[out] pk        - the public key.
/// @param[out] sk        - the secret key.
/// @param[in]  pk_seed   - seed for generating parts of public key.
/// @param[in]  sk_seed   - seed for generating secret key.
///
void generate_keypair_cyclic( cpk_t * pk, sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed );

///
/// @brief Generate compressed key pairs for cyclic rainbow.
///
/// @param[out] pk        - the public key.
/// @param[out] sk        - the compressed secret key.
/// @param[in]  pk_seed   - seed for generating parts of the public key.
/// @param[in]  sk_seed   - seed for generating the secret key.
///
void generate_compact_keypair_cyclic( cpk_t * pk, csk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed );

////////////////////////////////////

///
/// @brief Generate secret key for classic rainbow.
///
/// @param[out] sk        - the secret key.
/// @param[in]  sk_seed   - seed for generating the secret key.
///
void generate_secretkey(sk_t *sk, const unsigned char *sk_seed);

///
/// @brief Generate secret key for cyclic rainbow.
///
/// @param[out] sk        - the secret key.
/// @param[in]  pk_seed   - seed for generating parts of the pbulic key.
/// @param[in]  sk_seed   - seed for generating the secret key.
///
void generate_secretkey_cyclic( sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed );

////////////////////////////////////

///
/// @brief converting formats of public keys : from cyclic version to classic key
///
/// @param[out] pk       - the classic public key.
/// @param[in]  cpk      - the cyclic  public key.
///
void cpk_to_pk( pk_t * pk , const cpk_t * cpk );

void generate_identity_hash(unsigned char *digest, const unsigned char *id);

void multiply_identity_sk(sk_t *usk, const unsigned char *id_hash, sk_t *msk);

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
