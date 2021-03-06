///
///  @file api.h
///
///  Created by Bassham, Lawrence E (Fed) on 9/6/17.
///  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
///
///
///
///   @brief This is a sample 'api.h' for use 'sign.c'
///

#ifndef api_h
#define api_h


#include "rainbow_config.h"
#include "rainbow_keypair.h"

#include "api_config.h"


//  Set these three values apropriately for your algorithm


#if defined _RAINBOW_CLASSIC

#define CRYPTO_MASTER_SECRET_KEY_BYTES sizeof(msk_t)
#define CRYPTO_MASTER_PUBLIC_KEY_BYTES sizeof(mpk_t)

#define CRYPTO_USER_SECRET_KEY_BYTES sizeof(usk_t)
#define CRYPTO_USER_PUBLIC_KEY_BYTES sizeof(upk_t)

#elif defined _RAINBOW_CYCLIC

#define CRYPTO_MASTER_SECRET_KEY_BYTES sizeof(msk_t)
#define CRYPTO_MASTER_PUBLIC_KEY_BYTES sizeof(cpk_t)

#elif defined _RAINBOW_CYCLIC_COMPRESSED

#define CRYPTO_MASTER_SECRET_KEY_BYTES sizeof(csk_t)
#define CRYPTO_MASTER_PUBLIC_KEY_BYTES sizeof(cpk_t)

#else
error here
#endif


#define CRYPTO_BYTES _SIGNATURE_BYTE

// Change the algorithm name
#define CRYPTO_ALGNAME _S_NAME _SUFFIX


#ifdef  __cplusplus
extern  "C" {
#endif


int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk);

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk);

#ifdef  __cplusplus
}
#endif

#endif /* api_h */
