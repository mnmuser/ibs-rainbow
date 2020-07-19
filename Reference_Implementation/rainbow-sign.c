///  @file rainbow-sign.c
///  @brief A command-line tool for signing a file.
///

#include <stdio.h>
#include <stdint.h>

#include "rainbow_config.h"

#include "utils.h"

#include "api.h"



int main( int argc , char ** argv )
{
    printf("%s\n", CRYPTO_ALGNAME);

    printf("usk size: %lu\n", CRYPTO_USER_SECRET_KEY_BYTES);
    printf("upk size: %lu\n", CRYPTO_USER_PUBLIC_KEY_BYTES);
    printf("hash size: %d\n", _HASH_LEN);
    printf("signature size: %d\n\n", CRYPTO_BYTES);

    if (3 != argc) {
        printf("Usage:\n\n\trainbow-sign sk_file_name file_to_be_signed\n\n");
        return -1;
    }

    uint8_t *_sk = (uint8_t *) malloc(CRYPTO_USER_SECRET_KEY_BYTES);

    FILE *fp;
    int r = 0;

    fp = fopen(argv[1], "r");
    if (NULL == fp) {
        printf("fail to open secret key file.\n");
        return -1;
    }
    r = byte_fget(fp, _sk, CRYPTO_USER_SECRET_KEY_BYTES);
    fclose(fp);
    if (CRYPTO_USER_SECRET_KEY_BYTES != r) {
        printf("fail to load key file.\n");
        return -1;
    }

    unsigned char *msg = NULL;
    unsigned long long mlen = 0;
    r = byte_read_file(&msg, &mlen, argv[2]);
    if (0 != r) {
        printf("fail to read message file.\n");
        return -1;
	}

	unsigned char * signature = malloc( mlen + CRYPTO_BYTES );
	if( NULL == signature ) {
		printf("alloc memory for signature buffer fail.\n");
		return -1;
	}

	unsigned long long smlen = 0;
	r = crypto_sign( signature, &smlen, msg , mlen , _sk );
	if( 0 != r ) {
        printf("sign() fail.\n");

        free(msg);
        free(signature);
        free(_sk);

        return -1;
    }

	byte_fdump( stdout , CRYPTO_ALGNAME " signature"  , signature + mlen , CRYPTO_BYTES );
	printf("\n");

	free( msg );
	free( signature );
	free( _sk );

	return 0;
}

