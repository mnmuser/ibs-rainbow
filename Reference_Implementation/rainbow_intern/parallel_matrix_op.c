///  @file parallel_matrix_op.c
///  @brief the standard implementations for functions in parallel_matrix_op.h
///
///  the standard implementations for functions in parallel_matrix_op.h
///


#include <stdio.h>
#include "blas_comm.h"
#include "blas.h"

#include "parallel_matrix_op.h"
#include "rainbow_keypair.h"

#include "string.h"

#include "polynomial.h"
#include "rainbow_keypair_computation.h"


////////////////    Section: triangle matrix <-> rectangle matrix   ///////////////////////////////////


void UpperTrianglize(unsigned char *btriC, const unsigned char *bA, unsigned Awidth, unsigned size_batch) {
    unsigned char *runningC = btriC;
    unsigned Aheight = Awidth;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < i; j++) {
            unsigned idx = idx_of_trimat(j, i, Aheight);
            gf256v_add(btriC + idx * size_batch, bA + size_batch * (i * Awidth + j), size_batch);
        }
        gf256v_add(runningC, bA + size_batch * (i * Awidth + i), size_batch * (Aheight - i));
        runningC += size_batch * (Aheight - i);
    }
}


void quartic_UpperTrianglize(unsigned char *btriC, const unsigned char *bA, unsigned Awidth, unsigned size_batch) {

    unsigned char *runningC = btriC;
    unsigned Aheight = Awidth;

    unsigned char tmp_summand_A[(N_LINEAR_POLY + 1) / 2];
    unsigned char tmp_summand_B[(N_CUBIC_POLY + 1) / 2];
    unsigned char tmp_sum[N_QUARTIC_POLY];

    unsigned final_o = 0;

    unsigned final_e[15];

    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < i; j++) {
            unsigned idx = idx_of_trimat(j, i, Aheight);
            for (unsigned k = 0; k < size_batch * 2; k++) { //*2 because GF16

//                gf16_grade_n_poly_copy(tmp_summand_A, 0, btriC + idx * size_batch * N_QUARTIC_POLY,
//                                       N_QUARTIC_POLY * k, 1);

                gf16_grade_n_poly_copy(tmp_summand_B, 0, bA + size_batch * (i * Awidth + j) * N_QUARTIC_POLY,
                                       N_QUARTIC_POLY * k, 3);

                polynomial_add(btriC + idx * size_batch * N_QUARTIC_POLY, N_QUARTIC_POLY * k, 1, tmp_summand_B, 0,
                               N_CUBIC_POLY, _full_e_power2);

//                gf16_grade_n_poly_copy(btriC + idx * size_batch * N_QUARTIC_POLY, N_QUARTIC_POLY * k, tmp_sum,
//                                       0, 3);

                ///not working for layer 2:
                //gf16_quartic_poly_copy(btriC + idx * size_batch * N_QUARTIC_POLY(_ID), N_QUARTIC_POLY(_ID) * k,
                //                       bA + size_batch * (i * Awidth + j) * N_QUARTIC_POLY(_ID),
                //                       N_QUARTIC_POLY(_ID) * k);
            }
            //gf256v_add( btriC + idx*size_batch , bA + size_batch*(i*Awidth+j) , size_batch );
        }
        for (unsigned l = 0; l < size_batch * (Aheight - i) * 2; l++) {
//            gf16_grade_n_poly_copy(tmp_summand_A, 0, runningC, N_QUARTIC_POLY * l, 1);

            gf16_grade_n_poly_copy(tmp_summand_B, 0, bA + size_batch * (i * Awidth + i) * N_QUARTIC_POLY,
                                   l * N_QUARTIC_POLY, 3);

            polynomial_add(runningC, N_QUARTIC_POLY * l, 1, tmp_summand_B, 0, N_CUBIC_POLY, _full_e_power2);

//            gf16_grade_n_poly_copy(runningC, N_QUARTIC_POLY * l, tmp_sum, 0, 3);

//            gf16_quartic_poly_copy(runningC, N_QUARTIC_POLY(_ID) * l,
//                                   bA + size_batch * (i * Awidth + i) * N_QUARTIC_POLY(_ID),
//                                   l * N_QUARTIC_POLY(_ID));
        }
        //gf256v_add(runningC, bA + size_batch * (i * Awidth + i), size_batch * (Aheight - i)); /// ATTENTION GF256
        runningC += size_batch * (Aheight - i) * N_QUARTIC_POLY;
    }
}




/////////////////  Section: matrix multiplications  ///////////////////////////////

/// Q2, F1,  T1, _V1, _V1_BYTE, _O1, _O1_BYTE
void batch_trimat_madd_gf16(unsigned char *bC, const unsigned char *btriA,
                            const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth,
                            unsigned size_batch) {
    unsigned Awidth = Bheight;
    unsigned Aheight = Awidth;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (k < i) continue;
                gf16v_madd(bC, &btriA[(k - i) * size_batch], gf16v_get_ele(&B[j * size_Bcolvec], k), size_batch);
            }
            bC += size_batch;
        }
        btriA += (Aheight - i) * size_batch;
    }
}


/// F2 += F1 * T1
/// bC += BtriA * B
/// Bcolvec :: Größe eines Spaltenvektors in Matrix B
/// size_batch :: number of the batched elements in the corresponding position of the matrix.
/// -> _V1, _V1_BYTE, _O1, _O1_BYTE
void
quartic_batch_trimat_madd_gf16(unsigned char *bC, const unsigned char *btriA, const unsigned char *B, unsigned Bheight,
                               unsigned size_Bcolvec, unsigned Bwidth,
                               unsigned size_batch) {

    unsigned Awidth = Bheight;
    unsigned Aheight = Awidth;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (k < i) continue;
                quartic_gf16v_madd(bC, btriA, k - i, B, j, k, size_batch, size_Bcolvec);
            }
            bC += (size_batch * N_QUARTIC_POLY);
        }
        btriA += (Aheight - i) * _ID * size_batch;
    }
}


void batch_trimat_madd_gf256( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Awidth = Bheight;
    unsigned Aheight = Awidth;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(k<i) continue;
                gf256v_madd( bC , & btriA[ (k-i)*size_batch ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
        btriA += (Aheight-i)*size_batch;
    }
}

void quartic_batch_trimatTr_madd_gf16(unsigned char *bC, const unsigned char *btriA,
                                      const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth,
                                      unsigned size_batch) {
    unsigned Aheight = Bheight;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (i < k) continue;
                quartic_gf16v_madd2(bC, btriA, (idx_of_trimat(k, i, Aheight)), 1, B, j * size_Bcolvec, k,
                                    size_batch,
                                    size_Bcolvec); //TODO: sizeBcolvec?
                //gf16v_madd( bC , & btriA[ size_batch*(idx_of_trimat(k,i,Aheight)) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch * N_QUARTIC_POLY;
        }
    }
}

void batch_trimatTr_madd_gf16( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i<k) continue;
                gf16v_madd( bC , & btriA[ size_batch*(idx_of_trimat(k,i,Aheight)) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_trimatTr_madd_gf256( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i<k) continue;
                gf256v_madd( bC , & btriA[ size_batch*(idx_of_trimat(k,i,Aheight)) ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}




void batch_2trimat_madd_gf16( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i==k) continue;
                gf16v_madd( bC , & btriA[ size_batch*(idx_of_2trimat(i,k,Aheight)) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_2trimat_madd_gf256( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i==k) continue;
                gf256v_madd( bC , & btriA[ size_batch*(idx_of_2trimat(i,k,Aheight)) ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}



void batch_matTr_madd_gf16( unsigned char * bC , const unsigned char* A_to_tr , unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
        const unsigned char* bB, unsigned Bwidth, unsigned size_batch ) {
    unsigned Atr_height = Awidth;
    unsigned Atr_width = Aheight;
    for (unsigned i = 0; i < Atr_height; i++) {
        for (unsigned j = 0; j < Atr_width; j++) {
            gf16v_madd(bC, &bB[j * Bwidth * size_batch], gf16v_get_ele(&A_to_tr[size_Acolvec * i], j),
                       size_batch * Bwidth);
        }
        bC += size_batch * Bwidth;
    }
}

void
quartic_batch_matTr_madd_gf16(unsigned char *bC, const unsigned char *A_to_tr, unsigned Aheight, unsigned size_Acolvec,
                              unsigned Awidth,
                              const unsigned char *bB, unsigned Bwidth, unsigned size_batch) {
    unsigned Atr_height = Awidth;
    unsigned Atr_width = Aheight;
    for (unsigned i = 0; i < Atr_height; i++) {
        for (unsigned j = 0; j < Atr_width; j++) {
            quartic_gf16v_madd2(bC, bB, j, 0, A_to_tr, i, j, size_batch * Bwidth, size_Acolvec);
            //gf16v_madd(bC, &bB[j * Bwidth * size_batch], gf16v_get_ele(&A_to_tr[size_Acolvec * i], j),size_batch * Bwidth);
        }
        bC += size_batch * Bwidth * N_QUARTIC_POLY;
    }
}

void batch_matTr_madd_gf256(unsigned char *bC, const unsigned char *A_to_tr, unsigned Aheight, unsigned size_Acolvec,
                            unsigned Awidth,
                            const unsigned char *bB, unsigned Bwidth, unsigned size_batch) {
    unsigned Atr_height = Awidth;
    unsigned Atr_width = Aheight;
    for (unsigned i = 0; i < Atr_height; i++) {
        for (unsigned j = 0; j < Atr_width; j++) {
            gf256v_madd(bC, &bB[j * Bwidth * size_batch], gf256v_get_ele(&A_to_tr[size_Acolvec * i], j),
                        size_batch * Bwidth);
        }
        bC += size_batch*Bwidth;
    }
}


void quartic_batch_bmatTr_madd_gf16(unsigned char *bC, const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
                                    const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth,
                                    unsigned size_batch) {
    const unsigned char *bA = bA_to_tr;
    unsigned Aheight = Awidth_before_tr;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                quartic_gf16v_madd(bC, bA, i + k + Aheight, B, j, k, size_batch, size_Bcolvec);
                //gf16v_madd( bC , & bA[ size_batch*(i+k*Aheight) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch * N_QUARTIC_POLY;
        }
    }
}

void batch_bmatTr_madd_gf16( unsigned char *bC , const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch )
{
    const unsigned char *bA = bA_to_tr;
    unsigned Aheight = Awidth_before_tr;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf16v_madd( bC , & bA[ size_batch*(i+k*Aheight) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_bmatTr_madd_gf256( unsigned char *bC , const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch )
{
    const unsigned char *bA = bA_to_tr;
    unsigned Aheight = Awidth_before_tr;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf256v_madd( bC , & bA[ size_batch*(i+k*Aheight) ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}


void quartic_batch_mat_madd_gf16(unsigned char *bC, const unsigned char *bA, unsigned Aheight,
                                 const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth,
                                 unsigned size_batch) {
    unsigned Awidth = Bheight;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                quartic_gf16v_madd(bC, bA, k, B, j, k, size_batch, size_Bcolvec);
                //gf16v_madd( bC , & bA[ k*size_batch ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch * N_QUARTIC_POLY;
        }
        bA += (Awidth) * size_batch * _ID;
    }
}

void batch_mat_madd_gf16( unsigned char * bC , const unsigned char* bA , unsigned Aheight,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Awidth = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf16v_madd( bC , & bA[ k*size_batch ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
        bA += (Awidth)*size_batch;
    }
}

void batch_mat_madd_gf256( unsigned char * bC , const unsigned char* bA , unsigned Aheight,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Awidth = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf256v_madd( bC , & bA[ k*size_batch ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
        bA += (Awidth)*size_batch;
    }
}




////////////////////  Section: "quadratric" matrix evaluation  ///////////////////////////////



void batch_quad_trimat_eval_gf16( unsigned char * y, const unsigned char * trimat, const unsigned char * x, unsigned dim , unsigned size_batch )
{
///
///    assert( dim <= 128 );
///    assert( size_batch <= 128 );
    unsigned char tmp[256];

    unsigned char _x[256];
    for (unsigned i = 0; i < dim; i++) {
        _x[i] = gf16v_get_ele(x, i);
    }
    gf256v_set_zero(y, size_batch);
    for (unsigned i = 0; i < dim; i++) {
        gf256v_set_zero(tmp, size_batch);
        for (unsigned j = i; j < dim; j++) {
            gf16v_madd(tmp, trimat, _x[j], size_batch);
            trimat += size_batch;
        }
        gf16v_madd(y, tmp, _x[i], size_batch);
    }
}

void batch_quad_trimat_eval_gf256( unsigned char * y, const unsigned char * trimat, const unsigned char * x, unsigned dim , unsigned size_batch )
{
///
///    assert( dim <= 256 );
///    assert( size_batch <= 256 );
    unsigned char tmp[256];

    unsigned char _x[256];
    for(unsigned i=0;i<dim;i++) _x[i] = gf256v_get_ele( x , i );

    gf256v_set_zero( y , size_batch );
    for(unsigned i=0;i<dim;i++) {
        gf256v_set_zero( tmp , size_batch );
        for(unsigned j=i;j<dim;j++) {
           gf256v_madd( tmp , trimat , _x[j] , size_batch );
           trimat += size_batch;
        }
        gf256v_madd( y , tmp , _x[i] , size_batch );
    }
}







void batch_quad_recmat_eval_gf16( unsigned char * z, const unsigned char * y, unsigned dim_y, const unsigned char * mat,
        const unsigned char * x, unsigned dim_x , unsigned size_batch )
{
///
///    assert( dim_x <= 128 );
///    assert( dim_y <= 128 );
///    assert( size_batch <= 128 );
    unsigned char tmp[128];

    unsigned char _x[128];
    for(unsigned i=0;i<dim_x;i++) _x[i] = gf16v_get_ele( x , i );
    unsigned char _y[128];
    for(unsigned i=0;i<dim_y;i++) _y[i] = gf16v_get_ele( y , i );

    gf256v_set_zero( z , size_batch );
    for(unsigned i=0;i<dim_y;i++) {
        gf256v_set_zero( tmp , size_batch );
        for(unsigned j=0;j<dim_x;j++) {
           gf16v_madd( tmp , mat , _x[j] , size_batch );
           mat += size_batch;
        }
        gf16v_madd( z , tmp , _y[i] , size_batch );
    }
}


void batch_quad_recmat_eval_gf256( unsigned char * z, const unsigned char * y, unsigned dim_y, const unsigned char * mat,
        const unsigned char * x, unsigned dim_x , unsigned size_batch )
{
///
///    assert( dim_x <= 128 );
///    assert( dim_y <= 128 );
///    assert( size_batch <= 128 );
    unsigned char tmp[128];

    unsigned char _x[128];
    for(unsigned i=0;i<dim_x;i++) _x[i] = gf256v_get_ele( x , i );
    unsigned char _y[128];
    for(unsigned i=0;i<dim_y;i++) _y[i] = gf256v_get_ele( y , i );

    gf256v_set_zero( z , size_batch );
    for (unsigned i = 0; i < dim_y; i++) {
        gf256v_set_zero(tmp, size_batch);
        for (unsigned j = 0; j < dim_x; j++) {
            gf256v_madd(tmp, mat, _x[j], size_batch);
            mat += size_batch;
        }
        gf256v_madd(z, tmp, _y[i], size_batch);
    }
}

//TODO: rewrite universal quartic_gf16v_madd

void quartic_gf16v_madd(uint8_t *C, const uint8_t *A, unsigned A_pointer_index, const unsigned char *B,
                        unsigned B_pointer_index, unsigned B_offset, unsigned size_batch, unsigned size_Bcolvec) {

    ///SHOULD BE DONE BETTER (WIP)--///
    unsigned char tmp_product[(N_QUARTIC_POLY + 1) / 2]; // could be better calculated with i4.. in poly.c

    unsigned tmp_e[15]; //size is too big..

    unsigned tmp_o = 0;

    ///--SHOULD BE DONE BETTER (WIP)///

    for (unsigned l = 0; l < size_batch * 2; l++) { // *2 for gf16 (size is in byte)
        //the inner loop of gf16vmadd
        polynomial_mul(2, &A[(A_pointer_index) * _ID * size_batch], l * _ID, _lin_e_power2, 2,
                       &B[B_pointer_index * size_Bcolvec],
                       B_offset * _ID,
                       _lin_e_power2, &tmp_o, tmp_product, 0, tmp_e);


        polynomial_add(C, (l * N_QUARTIC_POLY), 1, tmp_product, 0, tmp_o, tmp_e);

    }
}


void quartic_gf16v_madd2(uint8_t *C, const uint8_t *Av, unsigned A_pointer_index, char A_linear, const unsigned char *B,
                         unsigned B_pointer_index, unsigned B_offset, unsigned size_batch,
                         unsigned size_Bcolvec) {

    ///SHOULD BE DONE BETTER (WIP)--///
    unsigned char tmp_product[(N_QUARTIC_POLY + 5) / 2]; // could be better calculated with i4.. in poly.c

    unsigned tmp_e[N_CUBIC_POLY + 2]; //size is too big..

    unsigned tmp_o = 0;

    unsigned const *e_A;
    unsigned o_A;
    unsigned A_loop_offset;

    unsigned C_grade;

    if (A_linear) {
        e_A = _lin_e_power2;
        o_A = N_LINEAR_POLY - 1;
        A_loop_offset = _ID;
        C_grade = 3;
    } else {
        e_A = _full_e_power2;
        o_A = N_QUADRATIC_POLY;
        A_loop_offset = N_QUARTIC_POLY;
        C_grade = 2;
    }
    ///--SHOULD BE DONE BETTER (WIP)///

    for (unsigned l = 0; l < size_batch * 2; l++) { // *2 for gf16 (size is in byte)
        //the inner loop of gf16vmadd

//        polynomial_print(o_A, &Av[(A_pointer_index) * size_batch * A_loop_offset], l * A_loop_offset, e_A, "l1_Q2:");
//        polynomial_print(2, &B[B_pointer_index * size_Bcolvec], B_offset * _ID, _lin_e_power2, "T1");

        polynomial_mul(o_A, &Av[(A_pointer_index) * size_batch * A_loop_offset], l * A_loop_offset, e_A, 2,
                       &B[B_pointer_index * size_Bcolvec],
                       B_offset * _ID,
                       _lin_e_power2, &tmp_o, tmp_product, 0, tmp_e);

//        polynomial_print(tmp_o, tmp_product, 0, tmp_e, "Produkt:");

//        polynomial_print(10,tmp_summand,(l * N_QUARTIC_POLY(_ID)),_full_e_power2,"tmp_sum");

        polynomial_add(C, (l * N_QUARTIC_POLY), C_grade, tmp_product, 0, tmp_o, tmp_e);

//        polynomial_print(15,C,(l * N_QUARTIC_POLY(_ID)),_full_e_power2,"Written:");
    }
}

// C is temp, A is S1, B is Qx
void quartic_gf16v_madd_to_grade(uint8_t *C, const uint8_t *A, unsigned A_pointer_index, const unsigned char *B,
                                 unsigned B_pointer_index, unsigned B_offset, unsigned B_grade, unsigned size_batch,
                                 unsigned size_Bcolvec) {

    ///SHOULD BE DONE BETTER (WIP)--///
    unsigned char tmp_product[(N_QUARTIC_POLY + 5) / 2];

    unsigned tmp_e[N_QUARTIC_POLY + 5]; //size is too big..

    unsigned tmp_o = 0;

    unsigned o2 = _grade_n_poly_terms(B_grade);

    ///--SHOULD BE DONE BETTER (WIP)///

    //A is grade 1, B is grade 3, C is empty

    for (unsigned l = 0; l < size_batch * 2; l++) { // *2 for gf16 (size is in byte)
        //the inner loop of gf16vmadd
        polynomial_mul(_ID, &A[(A_pointer_index) * _ID * size_batch], l * _ID, _lin_e_power2, o2,
                       &B[B_pointer_index * size_Bcolvec * N_QUARTIC_POLY],
                       B_offset * N_QUARTIC_POLY,
                       _full_e_power2, &tmp_o, tmp_product, 0, tmp_e);

//        polynomial_print(15,tmp_product,0,tmp_e,"Product:");

        polynomial_add(C, l * N_QUARTIC_POLY, 1, tmp_product, 0, tmp_o, tmp_e);

//        polynomial_print(15, C, (l * N_QUARTIC_POLY), final_e, "Sum:");

    }
}

//TODO: ONE gf16v_madd and you can choose grade of all inputs
//LESSONS LEARNED: don't try to unify to early

void calculate_values_public_key(unsigned char *upk, unsigned char *mpk, unsigned char *id) {
    unsigned char value_i;
    for (unsigned i = 0; i < sizeof(upk_t) * 2; i++) {
        value_i = polynomial_value(N_QUARTIC_POLY, mpk, i * N_QUARTIC_POLY, _full_e_power2, id);
        gf16v_set_ele(upk, i, value_i);
    }
}

void calculate_values_secret_key(unsigned char *usk, unsigned char *msk, unsigned char *id) {
    memcpy(usk, msk, LEN_SKSEED); //TODO: SEC: we don't want the random seed in the user-key
    usk += LEN_SKSEED;
    msk += LEN_SKSEED;

    for (unsigned i = 0; i < (sizeof(usk_t) - LEN_SKSEED) * 2; i++) {
        gf16v_set_ele(usk, i, polynomial_value(_ID, msk, i * _ID, _lin_e_power2, id));
    }
}


