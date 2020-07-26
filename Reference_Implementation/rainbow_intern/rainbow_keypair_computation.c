/// @file rainbow_keypair_computation.c
/// @brief Implementations for functions in rainbow_keypair_computation.h
///

#include "rainbow_keypair.h"
#include "rainbow_keypair_computation.h"

#include "blas_comm.h"
#include "rainbow_blas.h"
#include "polynomial.h"

#include <stdlib.h>
#include <string.h>


////////////////////////////////////////////////////////////////




void quartic_extcpk_to_pk(mpk_t *pk, const ext_mpk_t *cpk) {

    ///ATTENTION: look at those for-loops!

    const unsigned char *idx_l1 = cpk->l1_Q1;
    const unsigned char *idx_l2 = cpk->l2_Q1;
    for (unsigned i = 0; i < _V1; i++) {
        for (unsigned j = i; j < _V1; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY], idx_l1, _O1_BYTE * N_QUARTIC_POLY);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY]) + _O1_BYTE * N_QUARTIC_POLY, idx_l2,
                   _O2_BYTE * N_QUARTIC_POLY);
            idx_l1 += _O1_BYTE * N_QUARTIC_POLY;
            idx_l2 += _O2_BYTE * N_QUARTIC_POLY;
        }
    }

//    polynomial_print(15,pk->pk,0,_full_e_power2,"PK Q1:");
//    polynomial_print(15,pk->pk,15,_full_e_power2,"PK Q1:");
//    polynomial_print(15,pk->pk,30,_full_e_power2,"PK Q1:");
//    polynomial_print(15,pk->pk,45,_full_e_power2,"PK Q1:");
//    polynomial_print(15,pk->pk,60,_full_e_power2,"PK Q1:");
//    polynomial_print(15,pk->pk,126720 * 2 - 15,_full_e_power2,"PK Q1 (end):");

    idx_l1 = cpk->l1_Q2;
    idx_l2 = cpk->l2_Q2;
    for (unsigned i = 0; i < _V1; i++) {
        for (unsigned j = _V1; j < _V1 + _O1; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY], idx_l1, _O1_BYTE * N_QUARTIC_POLY);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY]) + _O1_BYTE * N_QUARTIC_POLY, idx_l2,
                   _O2_BYTE * N_QUARTIC_POLY);
            idx_l1 += _O1_BYTE * N_QUARTIC_POLY;
            idx_l2 += _O2_BYTE * N_QUARTIC_POLY;
        }
    }
    idx_l1 = cpk->l1_Q3;
    idx_l2 = cpk->l2_Q3;
    for (unsigned i = 0; i < _V1; i++) {
        for (unsigned j = _V1 + _O1; j < _PUB_N; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY], idx_l1, _O1_BYTE * N_QUARTIC_POLY);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY]) + _O1_BYTE * N_QUARTIC_POLY, idx_l2,
                   _O2_BYTE * N_QUARTIC_POLY);
            idx_l1 += _O1_BYTE * N_QUARTIC_POLY;
            idx_l2 += _O2_BYTE * N_QUARTIC_POLY;
        }
    }
    idx_l1 = cpk->l1_Q5;
    idx_l2 = cpk->l2_Q5;
    for (unsigned i = _V1; i < _V1 + _O1; i++) {
        for (unsigned j = i; j < _V1 + _O1; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY], idx_l1, _O1_BYTE * N_QUARTIC_POLY);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY]) + _O1_BYTE * N_QUARTIC_POLY, idx_l2,
                   _O2_BYTE * N_QUARTIC_POLY);
            idx_l1 += _O1_BYTE * N_QUARTIC_POLY;
            idx_l2 += _O2_BYTE * N_QUARTIC_POLY;
        }
    }
    idx_l1 = cpk->l1_Q6;
    idx_l2 = cpk->l2_Q6;
    for (unsigned i = _V1; i < _V1 + _O1; i++) {
        for (unsigned j = _V1 + _O1; j < _PUB_N; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY], idx_l1, _O1_BYTE * N_QUARTIC_POLY);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY]) + _O1_BYTE * N_QUARTIC_POLY, idx_l2,
                   _O2_BYTE * N_QUARTIC_POLY);
            idx_l1 += _O1_BYTE * N_QUARTIC_POLY;
            idx_l2 += _O2_BYTE * N_QUARTIC_POLY;
        }
    }
    idx_l1 = cpk->l1_Q9;
    idx_l2 = cpk->l2_Q9;
    for (unsigned i = _V1 + _O1; i < _PUB_N; i++) {
        for (unsigned j = i; j < _PUB_N; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY], idx_l1, _O1_BYTE * N_QUARTIC_POLY);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx * N_QUARTIC_POLY]) + _O1_BYTE * N_QUARTIC_POLY, idx_l2,
                   _O2_BYTE * N_QUARTIC_POLY);
            idx_l1 += _O1_BYTE * N_QUARTIC_POLY;
            idx_l2 += _O2_BYTE * N_QUARTIC_POLY;
        }
    }
}

void extcpk_to_pk(upk_t *pk, const ext_cpk_t *cpk) {
    const unsigned char *idx_l1 = cpk->l1_Q1;
    const unsigned char *idx_l2 = cpk->l2_Q1;
    for (unsigned i = 0; i < _V1; i++) {
        for (unsigned j = i; j < _V1; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx], idx_l1, _O1_BYTE);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx]) + _O1_BYTE, idx_l2, _O2_BYTE);
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q2;
    idx_l2 = cpk->l2_Q2;
    for (unsigned i = 0; i < _V1; i++) {
        for (unsigned j = _V1; j < _V1 + _O1; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx], idx_l1, _O1_BYTE);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx]) + _O1_BYTE, idx_l2, _O2_BYTE);
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q3;
    idx_l2 = cpk->l2_Q3;
    for (unsigned i = 0; i < _V1; i++) {
        for (unsigned j = _V1 + _O1; j < _PUB_N; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx], idx_l1, _O1_BYTE);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx]) + _O1_BYTE, idx_l2, _O2_BYTE);
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q5;
    idx_l2 = cpk->l2_Q5;
    for (unsigned i = _V1; i < _V1 + _O1; i++) {
        for (unsigned j = i; j < _V1 + _O1; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx], idx_l1, _O1_BYTE);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx]) + _O1_BYTE, idx_l2, _O2_BYTE);
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q6;
    idx_l2 = cpk->l2_Q6;
    for (unsigned i = _V1; i < _V1 + _O1; i++) {
        for (unsigned j = _V1 + _O1; j < _PUB_N; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx], idx_l1, _O1_BYTE);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx]) + _O1_BYTE, idx_l2, _O2_BYTE);
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q9;
    idx_l2 = cpk->l2_Q9;
    for (unsigned i = _V1 + _O1; i < _PUB_N; i++) {
        for (unsigned j = i; j < _PUB_N; j++) {
            unsigned pub_idx = idx_of_trimat(i, j, _PUB_N);
            memcpy(&pk->pk[_PUB_M_BYTE * pub_idx], idx_l1, _O1_BYTE);
            memcpy((&pk->pk[_PUB_M_BYTE * pub_idx]) + _O1_BYTE, idx_l2, _O2_BYTE);
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
}


/////////////////////////////////////////////////////////


static
void quartic_calculate_Q_from_F_ref(ext_mpk_t *cpk, const msk_t *sk) {
    //TODO: hey you, this could be done multi-threaded https://doi.org/10.1007/978-1-4419-5906-5 -> Related Work
/*
    Layer 1
    Computing :
    Q_pk.l1_F1s[i] = F_sk.l1_F1s[i]

    Q_pk.l1_F2s[i] = (F1* T1 + F2) + F1tr * t1
    Q_pk.l1_F5s[i] = UT( T1tr* (F1 * T1 + F2) )
*/
    const unsigned char *t2 = sk->t4;

    write_lin_to_quartic(cpk->l1_Q1, sk->l1_F1, _O1_BYTE * N_TRIANGLE_TERMS(_V1));

    write_lin_to_quartic(cpk->l1_Q2, sk->l1_F2, _O1_BYTE * _V1 * _O1);

    ///CHECK Q1
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q1, 0, _full_e_power2, "L1_Q1 ");
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q1, sizeof(cpk->l1_Q1) * 2 - N_QUARTIC_POLY, _full_e_power2, "L1_Q1 end");
    ///

    quartic_batch_trimat_madd(cpk->l1_Q2, sk->l1_F1, sk->t1, _V1, _V1_BYTE, _O1, _O1_BYTE); // Q2 += F1*T1

    set_quartic_zero(cpk->l1_Q3, _O1_BYTE * _V1 * _O2);
    set_quartic_zero(cpk->l1_Q5, _O1_BYTE * N_TRIANGLE_TERMS(_O1));
    set_quartic_zero(cpk->l1_Q6, _O1_BYTE * _O1 * _O2);
    set_quartic_zero(cpk->l1_Q9, _O1_BYTE * N_TRIANGLE_TERMS(_O2));

    // l1_Q5 : _O1_BYTE * _O1 * _O1
    // l1_Q9 : _O1_BYTE * _O2 * _O2
    // l2_Q5 : _O2_BYTE * _V1 * _O1
    // l2_Q9 : _O2_BYTE * _V1 * _O2
    unsigned size_tempQ = _O1_BYTE * _O1 * _O1 * N_QUARTIC_POLY;
    if (_O1_BYTE * _O2 * _O2 * N_QUARTIC_POLY > size_tempQ)
        size_tempQ = _O1_BYTE * _O2 * _O2 * N_QUARTIC_POLY;
    if (_O2_BYTE * _O1 * _O1 * N_QUARTIC_POLY > size_tempQ)
        size_tempQ = _O2_BYTE * _O1 * _O1 * N_QUARTIC_POLY;
    if (_O2_BYTE * _O2 * _O2 * N_QUARTIC_POLY > size_tempQ)
        size_tempQ = _O2_BYTE * _O2 * _O2 * N_QUARTIC_POLY;
    unsigned char *tempQ = (unsigned char *) aligned_alloc(32, size_tempQ + 32);

    memset(tempQ, 0, size_tempQ);  // l1_Q5


    quartic_batch_matTr_madd_gf16(tempQ, sk->t1, _V1, _V1_BYTE, _O1, cpk->l1_Q2, _O1, _O1_BYTE); // t1_tr*(F1*T1 + F2)

    polynomial_print(N_QUARTIC_POLY, tempQ, 0, _full_e_power2, "tempQ5(0): ");
    polynomial_print(N_QUARTIC_POLY, tempQ, size_tempQ * 2 - N_QUARTIC_POLY, _full_e_power2, "tempQ5(end): ");


    quartic_UpperTrianglize(cpk->l1_Q5, tempQ, _O1, _O1_BYTE);    // UT( ... )   // Q5 //TODO: copy?

    ///CHECK Q5
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q5, 0, _full_e_power2, "l1_Q5(0): ");
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q5, sizeof(cpk->l1_Q5) * 2 - N_QUARTIC_POLY, _full_e_power2, "L1_Q5 end");
    ///

    quartic_batch_trimatTr_madd_gf16(cpk->l1_Q2, sk->l1_F1, sk->t1, _V1, _V1_BYTE, _O1, _O1_BYTE); // Q2

    ///CHECK Q2
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q2, 0, _full_e_power2, "L1_Q2_1: ");
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q2, sizeof(cpk->l1_Q2) * 2 - N_QUARTIC_POLY, _full_e_power2,
                     "L1_Q2_2: "); //last position in Q2
    ///

/*
    Computing:
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2
    Q_pk.l1_F3s[i] =         F1_F1T_T2 + F2_T3                  ///linear
    Q_pk.l1_F6s[i] = T1tr* ( F1_F1T_T2 + F2_T3 ) + F2tr * t2    ///
    Q_pk.l1_F9s[i] = UT( T2tr* ( F1_T2 + F2_T3 ) )
*/
    quartic_batch_trimat_madd(cpk->l1_Q3, sk->l1_F1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE);         // F1*T2

    quartic_batch_mat_madd_gf16(cpk->l1_Q3, sk->l1_F2, _V1, sk->t3, _O1, _O1_BYTE, _O2, _O1_BYTE);   // F1_T2 + F2_T3

    set_quartic_zero(tempQ, _O1_BYTE * _O2 * _O2);                                              // l1_Q9

//    polynomial_print(15, tempQ, 0, _full_e_power2, "tmpQ9(zero): ");
//    polynomial_print(15, cpk->l1_Q3, 0, _full_e_power2, "tmp_Q3(0): ");
//    polynomial_print(_ID,t2,0,_lin_e_power2,"t2:");

    quartic_batch_matTr_madd_gf16(tempQ, t2, _V1, _V1_BYTE, _O2, cpk->l1_Q3, _O2,
                                  _O1_BYTE);           // T2tr * ( F1_T2 + F2_T3 )

    quartic_UpperTrianglize(cpk->l1_Q9, tempQ, _O2, _O1_BYTE);                                   // Q9

    ///CHECK Q9
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q9, 0, _full_e_power2, "l1_Q9(0): ");
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q9, sizeof(cpk->l1_Q9) * 2 - N_QUARTIC_POLY, _full_e_power2, "L1_Q9 end");
    ///

    quartic_batch_trimatTr_madd_gf16(cpk->l1_Q3, sk->l1_F1, t2, _V1, _V1_BYTE, _O2,
                                     _O1_BYTE);        // F1_F1T_T2 + F2_T3  // Q3

    ///CHECK Q3
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q3, 0, _full_e_power2, "Q3:");
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q3, sizeof(cpk->l1_Q3) * 2 - N_QUARTIC_POLY, _full_e_power2,
                     "L1_Q3 end: ");
    ///

    quartic_batch_bmatTr_madd_gf16(cpk->l1_Q6, sk->l1_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE);       // F2tr*T2

    quartic_batch_matTr_madd_gf16(cpk->l1_Q6, sk->t1, _V1, _V1_BYTE, _O1, cpk->l1_Q3, _O2,
                                  _O1_BYTE);    // Q6

    ///CHECK Q6
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q6, 0, _full_e_power2, "Q6:");
    polynomial_print(N_QUARTIC_POLY, cpk->l1_Q6, sizeof(cpk->l1_Q6) * 2 - N_QUARTIC_POLY, _full_e_power2,
                     "L1_Q6 end: ");
    ///


    /////////////////// LAYER 2 /////////////////////////

    /*
    layer 2
    Computing:
    Q1 = F1
    Q2 = F1_F1T*T1 + F2
    Q5 = UT( T1tr( F1*T1 + F2 )  + F5 )
*/
    write_lin_to_quartic(cpk->l2_Q1, sk->l2_F1, _O2_BYTE * N_TRIANGLE_TERMS(_V1));

    ///CHECK Q1
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q1, 0, _full_e_power2, "L2_Q1 ");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q1, sizeof(cpk->l2_Q1) * 2 - N_QUARTIC_POLY, _full_e_power2, "L2_Q1 end");
    ///

    write_lin_to_quartic(cpk->l2_Q2, sk->l2_F2, _O2_BYTE * _V1 * _O1);
    quartic_batch_trimat_madd(cpk->l2_Q2, sk->l2_F1, sk->t1, _V1, _V1_BYTE, _O1, _O2_BYTE);      // F1*T1 + F2

    write_lin_to_quartic(cpk->l2_Q5, sk->l2_F5, _O2_BYTE * N_TRIANGLE_TERMS(_O1));
    set_quartic_zero(tempQ, _O2_BYTE * _O1 * _O1);                                               // l2_Q5
    quartic_batch_matTr_madd_gf16(tempQ, sk->t1, _V1, _V1_BYTE, _O1, cpk->l2_Q2, _O1,
                                  _O2_BYTE);        // t1_tr*(F1*T1 + F2)
    quartic_UpperTrianglize(cpk->l2_Q5, tempQ, _O1,
                            _O2_BYTE);                             // UT( ... )   // Q5

    ///CHECK Q5
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q5, 0, _full_e_power2, "l2_Q5(0): ");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q5, sizeof(cpk->l2_Q5) * 2 - N_QUARTIC_POLY, _full_e_power2, "L2_Q5 end");
    ///

    quartic_batch_trimatTr_madd_gf16(cpk->l2_Q2, sk->l2_F1, sk->t1, _V1, _V1_BYTE, _O1, _O2_BYTE);    // Q2

    ///CHECK Q2
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q2, 0, _full_e_power2, "L2_Q2_1: ");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q2, sizeof(cpk->l2_Q2) * 2 - N_QUARTIC_POLY, _full_e_power2,
                     "L2_Q2_2: "); //last position in Q2
    ///

/*
    Computing:
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2

    Q3 =        F1_F1T*T2 + F2*T3 + F3
    Q9 = UT( T2tr*( F1*T2 + F2*T3 + F3 )  +      T3tr*( F5*T3 + F6 ) )
    Q6 = T1tr*( F1_F1T*T2 + F2*T3 + F3 )  + F2Tr*T2 + F5_F5T*T3 + F6
*/
    write_lin_to_quartic(cpk->l2_Q3, sk->l2_F3, _O2_BYTE * _V1 * _O2);
    quartic_batch_trimat_madd(cpk->l2_Q3, sk->l2_F1, t2, _V1, _V1_BYTE, _O2, _O2_BYTE);         // F1*T2 + F3
    quartic_batch_mat_madd_gf16(cpk->l2_Q3, sk->l2_F2, _V1, sk->t3, _O1, _O1_BYTE, _O2,
                                _O2_BYTE);   // F1_T2 + F2_T3 + F3

    set_quartic_zero(tempQ, _O2_BYTE * _O2 * _O2);                                              // l2_Q9
    quartic_batch_matTr_madd_gf16(tempQ, t2, _V1, _V1_BYTE, _O2, cpk->l2_Q3, _O2,
                                  _O2_BYTE);           // T2tr * ( ..... )

    write_lin_to_quartic(cpk->l2_Q6, sk->l2_F6, _O2_BYTE * _O1 * _O2);

    quartic_batch_trimat_madd(cpk->l2_Q6, sk->l2_F5, sk->t3, _O1, _O1_BYTE, _O2, _O2_BYTE);      // F5*T3 + F6
    quartic_batch_matTr_madd_gf16(tempQ, sk->t3, _O1, _O1_BYTE, _O2, cpk->l2_Q6, _O2,
                                  _O2_BYTE);       // T2tr*( ..... ) + T3tr*( ..... )

    set_quartic_zero(cpk->l2_Q9, _O2_BYTE * N_TRIANGLE_TERMS(_O2));
    quartic_UpperTrianglize(cpk->l2_Q9, tempQ, _O2, _O2_BYTE);                                   // Q9

    ///CHECK Q9
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q9, 0, _full_e_power2, "l2_Q9(0): ");

    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q9, sizeof(cpk->l2_Q9) * 2 - N_QUARTIC_POLY * 4, _full_e_power2,
                     "L2_Q9 end -45");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q9, sizeof(cpk->l2_Q9) * 2 - N_QUARTIC_POLY * 3, _full_e_power2,
                     "L2_Q9 end -30");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q9, sizeof(cpk->l2_Q9) * 2 - N_QUARTIC_POLY * 2, _full_e_power2,
                     "L2_Q9 end -15");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q9, sizeof(cpk->l2_Q9) * 2 - N_QUARTIC_POLY, _full_e_power2, "L2_Q9 end");
    ///

    quartic_batch_trimatTr_madd_gf16(cpk->l2_Q3, sk->l2_F1, t2, _V1, _V1_BYTE, _O2,
                                     _O2_BYTE);        // F1_F1T_T2 + F2_T3 + F3 // Q3

    ///CHECK Q3
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q3, 0, _full_e_power2, "L2_Q3:");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q3, sizeof(cpk->l2_Q3) * 2 - N_QUARTIC_POLY, _full_e_power2,
                     "L2_Q3 end: ");
    ///

    quartic_batch_bmatTr_madd_gf16(cpk->l2_Q6, sk->l2_F2, _O1, t2, _V1, _V1_BYTE, _O2,
                                   _O2_BYTE);       //  F5*T3 + F6 +  F2tr*T2
    quartic_batch_trimatTr_madd_gf16(cpk->l2_Q6, sk->l2_F5, sk->t3, _O1, _O1_BYTE, _O2,
                                     _O2_BYTE);    //   F2tr*T2 + F5_F5T*T3 + F6
    quartic_batch_matTr_madd_gf16(cpk->l2_Q6, sk->t1, _V1, _V1_BYTE, _O1, cpk->l2_Q3, _O2, _O2_BYTE);    // Q6

    ///CHECK Q6
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q6, 0, _full_e_power2, "l2_Q6:");
    polynomial_print(N_QUARTIC_POLY, cpk->l2_Q6, sizeof(cpk->l2_Q6) * 2 - N_QUARTIC_POLY, _full_e_power2,
                     "L2_Q6 end: ");
    ///

    //TODO: write polynomial_fill_checker?
    //TODO: GF16-naming
    memset(tempQ, 0, size_tempQ + 32);
    free(tempQ);
}

static
void calculate_Q_from_F_ref(ext_cpk_t *Qs, const usk_t *Fs, const usk_t *Ts) {
/*
    Layer 1
    Computing :
    Q_pk.l1_F1s[i] = F_sk.l1_F1s[i]

    Q_pk.l1_F2s[i] = (F1* T1 + F2) + F1tr * t1
    Q_pk.l1_F5s[i] = UT( T1tr* (F1 * T1 + F2) )
*/
    const unsigned char *t2 = Ts->t4;

    memcpy(Qs->l1_Q1, Fs->l1_F1, _O1_BYTE * N_TRIANGLE_TERMS(_V1));

    memcpy(Qs->l1_Q2, Fs->l1_F2, _O1_BYTE * _V1 * _O1);
    batch_trimat_madd(Qs->l1_Q2, Fs->l1_F1, Ts->t1, _V1, _V1_BYTE, _O1, _O1_BYTE);    // F1*T1 + F2

    memset(Qs->l1_Q3, 0, _O1_BYTE * _V1 * _O2);
    memset(Qs->l1_Q5, 0, _O1_BYTE * N_TRIANGLE_TERMS(_O1));
    memset(Qs->l1_Q6, 0, _O1_BYTE * _O1 * _O2);
    memset(Qs->l1_Q9, 0, _O1_BYTE * N_TRIANGLE_TERMS(_O2));

    // l1_Q5 : _O1_BYTE * _O1 * _O1
    // l1_Q9 : _O1_BYTE * _O2 * _O2
    // l2_Q5 : _O2_BYTE * _V1 * _O1
    // l2_Q9 : _O2_BYTE * _V1 * _O2
    unsigned size_tempQ = _O1_BYTE * _O1 * _O1;
    if (_O1_BYTE * _O2 * _O2 > size_tempQ) size_tempQ = _O1_BYTE * _O2 * _O2;
    if (_O2_BYTE * _O1 * _O1 > size_tempQ) size_tempQ = _O2_BYTE * _O1 * _O1;
    if (_O2_BYTE * _O2 * _O2 > size_tempQ) size_tempQ = _O2_BYTE * _O2 * _O2;
    unsigned char *tempQ = (unsigned char *) aligned_alloc(32, size_tempQ + 32);

    memset(tempQ, 0, _O1_BYTE * _O1 * _O1);   // l1_Q5
    batch_matTr_madd(tempQ, Ts->t1, _V1, _V1_BYTE, _O1, Qs->l1_Q2, _O1, _O1_BYTE);  // t1_tr*(F1*T1 + F2)
    UpperTrianglize(Qs->l1_Q5, tempQ, _O1, _O1_BYTE);    // UT( ... )   // Q5

    batch_trimatTr_madd(Qs->l1_Q2, Fs->l1_F1, Ts->t1, _V1, _V1_BYTE, _O1, _O1_BYTE);    // Q2
/*
    Computing:
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2
    Q_pk.l1_F3s[i] =         F1_F1T_T2 + F2_T3
    Q_pk.l1_F6s[i] = T1tr* ( F1_F1T_T2 + F2_T3 ) + F2tr * t2
    Q_pk.l1_F9s[i] = UT( T2tr* ( F1_T2 + F2_T3 ) )
*/
    batch_trimat_madd(Qs->l1_Q3, Fs->l1_F1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE);         // F1*T2
    batch_mat_madd(Qs->l1_Q3, Fs->l1_F2, _V1, Ts->t3, _O1, _O1_BYTE, _O2, _O1_BYTE);   // F1_T2 + F2_T3

    memset(tempQ, 0, _O1_BYTE * _O2 * _O2);                                              // l1_Q9
    batch_matTr_madd(tempQ, t2, _V1, _V1_BYTE, _O2, Qs->l1_Q3, _O2, _O1_BYTE);           // T2tr * ( F1_T2 + F2_T3 )
    UpperTrianglize(Qs->l1_Q9, tempQ, _O2, _O1_BYTE);                                   // Q9

    batch_trimatTr_madd(Qs->l1_Q3, Fs->l1_F1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE);        // F1_F1T_T2 + F2_T3  // Q3

    batch_bmatTr_madd(Qs->l1_Q6, Fs->l1_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE);       // F2tr*T2
    batch_matTr_madd(Qs->l1_Q6, Ts->t1, _V1, _V1_BYTE, _O1, Qs->l1_Q3, _O2, _O1_BYTE);    // Q6

/*
    layer 2
    Computing:
    Q1 = F1
    Q2 = F1_F1T*T1 + F2
    Q5 = UT( T1tr( F1*T1 + F2 )  + F5 )
*/
    memcpy(Qs->l2_Q1, Fs->l2_F1, _O2_BYTE * N_TRIANGLE_TERMS(_V1));

    memcpy(Qs->l2_Q2, Fs->l2_F2, _O2_BYTE * _V1 * _O1);
    batch_trimat_madd(Qs->l2_Q2, Fs->l2_F1, Ts->t1, _V1, _V1_BYTE, _O1, _O2_BYTE);      // F1*T1 + F2

    memcpy(Qs->l2_Q5, Fs->l2_F5, _O2_BYTE * N_TRIANGLE_TERMS(_O1));
    memset(tempQ, 0, _O2_BYTE * _O1 * _O1);                                               // l2_Q5
    batch_matTr_madd(tempQ, Ts->t1, _V1, _V1_BYTE, _O1, Qs->l2_Q2, _O1, _O2_BYTE);        // t1_tr*(F1*T1 + F2)
    UpperTrianglize(Qs->l2_Q5, tempQ, _O1, _O2_BYTE);                                     // UT( ... )   // Q5

    batch_trimatTr_madd(Qs->l2_Q2, Fs->l2_F1, Ts->t1, _V1, _V1_BYTE, _O1, _O2_BYTE);    // Q2

/*
    Computing:
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2

    Q3 =        F1_F1T*T2 + F2*T3 + F3
    Q9 = UT( T2tr*( F1*T2 + F2*T3 + F3 )  +      T3tr*( F5*T3 + F6 ) )
    Q6 = T1tr*( F1_F1T*T2 + F2*T3 + F3 )  + F2Tr*T2 + F5_F5T*T3 + F6
*/
    memcpy(Qs->l2_Q3, Fs->l2_F3, _O2_BYTE * _V1 * _O2);
    batch_trimat_madd(Qs->l2_Q3, Fs->l2_F1, t2, _V1, _V1_BYTE, _O2, _O2_BYTE);         // F1*T2 + F3
    batch_mat_madd(Qs->l2_Q3, Fs->l2_F2, _V1, Ts->t3, _O1, _O1_BYTE, _O2, _O2_BYTE);   // F1_T2 + F2_T3 + F3

    memset(tempQ, 0, _O2_BYTE * _O2 * _O2);                                              // l2_Q9
    batch_matTr_madd(tempQ, t2, _V1, _V1_BYTE, _O2, Qs->l2_Q3, _O2, _O2_BYTE);           // T2tr * ( ..... )

    memcpy(Qs->l2_Q6, Fs->l2_F6, _O2_BYTE * _O1 * _O2);

    batch_trimat_madd(Qs->l2_Q6, Fs->l2_F5, Ts->t3, _O1, _O1_BYTE, _O2, _O2_BYTE);      // F5*T3 + F6
    batch_matTr_madd(tempQ, Ts->t3, _O1, _O1_BYTE, _O2, Qs->l2_Q6, _O2,
                     _O2_BYTE);       // T2tr*( ..... ) + T3tr*( ..... )
    memset(Qs->l2_Q9, 0, _O2_BYTE * N_TRIANGLE_TERMS(_O2));
    UpperTrianglize(Qs->l2_Q9, tempQ, _O2, _O2_BYTE);                                   // Q9

    batch_trimatTr_madd(Qs->l2_Q3, Fs->l2_F1, t2, _V1, _V1_BYTE, _O2, _O2_BYTE);        // F1_F1T_T2 + F2_T3 + F3 // Q3

    batch_bmatTr_madd(Qs->l2_Q6, Fs->l2_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O2_BYTE);       //  F5*T3 + F6 +  F2tr*T2
    batch_trimatTr_madd(Qs->l2_Q6, Fs->l2_F5, Ts->t3, _O1, _O1_BYTE, _O2, _O2_BYTE);    //   F2tr*T2 + F5_F5T*T3 + F6
    batch_matTr_madd(Qs->l2_Q6, Ts->t1, _V1, _V1_BYTE, _O1, Qs->l2_Q3, _O2, _O2_BYTE);    // Q6

    memset(tempQ, 0, size_tempQ + 32);
    free(tempQ);
}

/////////////////////////////////////////////////////////////////


// Choosing implementations depends on the macros: _BLAS_SSE_ and _BLAS_AVX2_
#if defined(_BLAS_SSE_) || defined(_BLAS_AVX2_)
#include "rainbow_keypair_computation_simd.h"
#define quartic_calculate_Q_from_F_impl        calculate_Q_from_F_simd
#define calculate_F_from_Q_impl        calculate_F_from_Q_simd
#define calculate_Q_from_F_cyclic_impl calculate_Q_from_F_cyclic_simd
#else
#define quartic_calculate_Q_from_F_impl        quartic_calculate_Q_from_F_ref
#define calculate_Q_from_F_impl        calculate_Q_from_F_ref
#define calculate_F_from_Q_impl        calculate_F_from_Q_ref
#endif

void calculate_Q_from_F(ext_cpk_t *Qs, const usk_t *Fs, const usk_t *Ts) {
    calculate_Q_from_F_impl(Qs, Fs, Ts);
}

void quartic_calculate_Q_from_F(ext_mpk_t *Qs, const msk_t *Fs) {
    quartic_calculate_Q_from_F_impl(Qs, Fs);
}

void write_lin_to_quartic(unsigned char *q, const unsigned char *f, unsigned long length_f) {
    set_quartic_zero(q, length_f); // to have a clean polynom
    for (unsigned x = 0; x < length_f * 2; x++) {
        for (unsigned i = 0; i < N_LINEAR_POLY; i++) {
            gf16v_set_ele(q, i + (N_QUARTIC_POLY * x),
                          gf16v_get_ele(f, i + x * N_LINEAR_POLY)); // plus one, to write after the constant field
        }
    }
}

void write_lin_wo_const_to_quadratic(unsigned char *dest, const unsigned char *src, unsigned long length_f) {
    set_quartic_zero(dest, length_f); // to have a clean polynom
    for (unsigned x = 0; x < length_f * 2; x++) {
        for (unsigned i = 0; i < _ID; i++) {
            gf16v_set_ele(dest, i + (N_QUADRATIC_POLY * x),
                          gf16v_get_ele(src, i + x * _ID)); // plus one, to write after the constant field
        }
    }
}

void set_quartic_zero(unsigned char *q, const unsigned length) {
    memset(q, 0, length * N_QUARTIC_POLY);
}

/// copy a gf16-polynomial of grade n
/// \param dest
/// \param gf16_offset_dest
/// \param src
/// \param gf16_offset_src
/// \param grade_n
void gf16_grade_n_poly_copy(unsigned char *dest, unsigned gf16_offset_dest, const unsigned char *src,
                            unsigned gf16_offset_src, unsigned grade_n) {
    for (unsigned i = 0; i < _grade_n_poly_terms(grade_n); i++) {
        gf16v_set_ele(dest, i + gf16_offset_dest, gf16v_get_ele(src, gf16_offset_src + i));
    }
}

void gf16_copy(unsigned char *dest, unsigned gf16_offset_dest, const unsigned char *src,
               unsigned gf16_offset_src, unsigned size_values) {
    for (unsigned i = 0; i < size_values; i++) {
        gf16v_set_ele(dest, i + gf16_offset_dest, gf16v_get_ele(src, gf16_offset_src + i));
    }
}

