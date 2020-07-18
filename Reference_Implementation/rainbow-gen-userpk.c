///  @file rainbow-gen-userpk.c
///  @brief A command-line tool for generating the user public keys.

#include <rainbow_keypair.h>
#include <stdlib.h>
#include <stdio.h>
#include <api.h>
#include <utils.h>
#include <string.h>
#include <blas.h>

int main(int argc, char **argv) {

    printf("%s\n", CRYPTO_ALGNAME);

    printf("msk size: %lu\n", CRYPTO_MASTER_SECRET_KEY_BYTES);
    printf("mpk size: %lu\n", CRYPTO_MASTER_PUBLIC_KEY_BYTES);
    printf("hash size: %d\n", _HASH_LEN);
    printf("signature size: %d\n\n", CRYPTO_BYTES);

    if (4 != argc) {
        printf("Usage:\n\n\trainbow-gen-userpk mpk_file_name identity upk_file_name\n\n");
        return -1;
    }


    //malloc upk

    uint8_t *_mpk = malloc(CRYPTO_MASTER_PUBLIC_KEY_BYTES);
    FILE *fp;
    int r = 0; // TODO: better unsigned?

    fp = fopen( argv[1] , "r");
    if( NULL == fp ) {
        printf("fail to open master key file.\n");
        return -1;
    }
    r = byte_fget(fp, _mpk, CRYPTO_MASTER_PUBLIC_KEY_BYTES);
    fclose(fp);
    if (CRYPTO_MASTER_PUBLIC_KEY_BYTES != r) {
        printf("fail to load key file.\n");
        return -1;
    }

    unsigned char * _identity = malloc((_ID+1)/2);
    r = 0;

    generate_identity_hash(_identity,argv[2],strlen(argv[2]));
    if( NULL == _identity ) {
        printf("fail to create identity hash.\n");
        return -1;
    }
    printf("identity hash: ");
    for(unsigned i = 0; i<_ID; i++){
        printf("%hhu ",gf16v_get_ele(_identity,i));
    }
    printf("\n\n");

    uint8_t *_upk = malloc(sizeof(upk_t));

    //calculate usk with mpk and ID

    int re = calculate_upk(_upk,_mpk,_identity);
    if (0 != re) {
        printf("%s generate user-public-key fails.\n", CRYPTO_ALGNAME);
        return -1;
    }

    //write upk to disk
    fp = fopen(argv[3], "w+");
    if (NULL == fp) {
        printf("fail to open user public key file.\n");
        return -1;
    }
    byte_fdump(fp, CRYPTO_ALGNAME " user public key", _upk, sizeof(upk_t)); //pk speichern und beschriften
    fclose(fp);

    //free upk and rest

    free(_mpk);
    free(_identity);
    free(_upk);

    return 0;
}