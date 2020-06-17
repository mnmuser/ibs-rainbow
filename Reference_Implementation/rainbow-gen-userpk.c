///  @file rainbow-gen-userpk.c
///  @brief A command-line tool for generating the user public keys.

#include <rainbow_keypair.h>
#include <stdlib.h>
#include <stdio.h>
#include <api.h>

int main(int argc, char **argv) {

    printf("%s\n", CRYPTO_ALGNAME);

    printf("msk size: %lu\n", CRYPTO_SECRETKEYBYTES);
    printf("mpk size: %lu\n", CRYPTO_PUBLICKEYBYTES);
    printf("hash size: %d\n", _HASH_LEN);
    printf("signature size: %d\n\n", CRYPTO_BYTES);

    if (3 != argc) {
        printf("Usage:\n\n\trainbow-gen-userpk mpk_file_name identity\n\n");
        return -1;
    }


    //malloc upk

    uint8_t *_upk = malloc(sizeof(upk_t));
    FILE *fp;

    //calculate usk with mpk and ID

    ///calculate_upk

    //write upk to disk

    //free upk

    free(_upk);

    return 0;
}