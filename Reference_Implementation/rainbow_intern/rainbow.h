/// @file rainbow.h
/// @brief APIs for rainbow.
///

#ifndef _RAINBOW_H_
#define _RAINBOW_H_

#include "rainbow_config.h"
#include "rainbow_keypair.h"

#include <stdint.h>

#ifdef  __cplusplus
extern  "C" {
#endif


///
/// @brief Signing function for classical secret key.
///
/// @param[out] signature - the signature.
/// @param[in]  sk        - the secret key.
/// @param[in]  digest    - the digest.
///
int rainbow_sign(uint8_t *signature, const usk_t *sk, const uint8_t *digest);

///
/// @brief Verifying function.
///
/// @param[in]  digest    - the digest.
/// @param[in]  signature - the signature.
/// @param[in]  pk        - the public key.
/// @return 0 for successful verified. -1 for failed verification.
///
int rainbow_verify(const uint8_t *digest, const uint8_t *signature, const upk_t *pk);



#ifdef  __cplusplus
}
#endif


#endif // _RAINBOW_H_
