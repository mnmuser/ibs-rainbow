///  @file rainbow-gen-usersk.c
///  @brief A command-line tool for generating the user secret keys.

#include <rainbow_keypair.h>
#include <stdlib.h>

int main(int argc, char **argv) {

    //malloc usk

    usk_t *usk = malloc(sizeof(usk_t));

    //calculate usk with msk and ID

    ///calculate_usk

    //write usk to disk

    //free usk

    free(usk);

    return 0;
}