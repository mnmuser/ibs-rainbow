///  @file rainbow-gen-userpk.c
///  @brief A command-line tool for generating the user public keys.

#include <rainbow_keypair.h>
#include <stdlib.h>

int main(int argc, char **argv) {

    //malloc upk

    upk_t *upk = malloc(sizeof(upk_t));

    //calculate usk with mpk and ID

    ///calculate_upk

    //write upk to disk

    //free upk

    free(upk);

    return 0;
}