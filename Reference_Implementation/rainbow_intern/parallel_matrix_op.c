///  @file parallel_matrix_op.c
///  @brief the standard implementations for functions in parallel_matrix_op.h
///
///  the standard implementations for functions in parallel_matrix_op.h
///


#include "blas_comm.h"
#include "blas.h"

#include "parallel_matrix_op.h"
#include "rainbow_keypair.h"

#include "stdlib.h"

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

// A is grade 3, C is empty, should be modified A afterwards
void quartic_UpperTrianglize(unsigned char *btriC, const unsigned char *bA, unsigned Awidth, unsigned size_batch) {
    unsigned char *runningC = btriC;
    unsigned Aheight = Awidth;
    unsigned o = N_CUBIC_POLY(_ID);
    int e;
    int final_o;
    unsigned char *tmp_C = malloc(N_QUARTIC_POLY(_ID));
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < i; j++) {
            unsigned idx = idx_of_trimat(j, i, Aheight);
            for (unsigned k = 0; k < size_batch * 2; k++) { //gf256 operates on bytes
                //TODO: memcpy C, then go on ;)
                gf16_quartic_poly_copy(tmp_C, btriC, (2 * idx * size_batch) + k * N_QUARTIC_POLY(_ID));
                polynomial_add(N_CUBIC_POLY(_ID), tmp_C, full_e_power2, N_CUBIC_POLY(_ID), bA, full_e_power2, &final_o,
                               btriC, k * N_QUARTIC_POLY(_ID), &e);
            }
            //gf256v_add( btriC + idx*size_batch , bA + size_batch*(i*Awidth+j) , size_batch );
        }
        gf256v_add(runningC, bA + size_batch * (i * Awidth + i), size_batch * (Aheight - i));
        runningC += size_batch * (Aheight - i);
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
            bC += (size_batch * N_QUARTIC_POLY(_ID));
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
                quartic_gf16v_madd2(bC, btriA, (idx_of_trimat(k, i, Aheight)), 1, B, j * size_Bcolvec, k, size_batch,
                                    size_Bcolvec);
                //gf16v_madd( bC , & btriA[ size_batch*(idx_of_trimat(k,i,Aheight)) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch * N_QUARTIC_POLY(_ID);
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
        bC += size_batch * Bwidth * N_QUARTIC_POLY(_ID);
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
            bC += size_batch;
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
            bC += size_batch * N_QUARTIC_POLY(_ID);
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

void quartic_gf16v_madd(uint8_t *C, const uint8_t *A, unsigned A_pointer_index, const unsigned char *B,
                        unsigned B_pointer_index, unsigned B_offset, unsigned size_batch, unsigned size_Bcolvec) {

    ///SHOULD BE DONE BETTER (WIP)--///
    int e_linear[2] = {2, 3}; // the structure of the sk-fields (do we need a constant factor?)
    unsigned char tmp_product[(N_QUARTIC_POLY(_ID) + 1) / 2]; // could be better calculated with i4.. in poly.c
    unsigned char tmp_summand[(_ID + 2) / 2]; //GF16, round up, one extra field for constant

    int tmp_e[15]; //size is too big..
    int final_e[15];

    int tmp_o = 0;
    int final_o = 0;

    //#define e_standard(n) (n+1) -> mal ersetzen
    int full_e_power2[15] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    ///--SHOULD BE DONE BETTER (WIP)///

    for (unsigned l = 0; l < size_batch * 2; l++) { // *2 for gf16 (size is in byte)
        //the inner loop of gf16vmadd
        polynomial_mul(2, &A[(A_pointer_index) * _ID * size_batch], l, e_linear, 2, &B[B_pointer_index * size_Bcolvec],
                       B_offset * _ID, // u.U l * _ID ????
                       e_linear, &tmp_o, tmp_product, 0, tmp_e);

        gf16_lin_poly_copy(tmp_summand, C, (l * N_QUARTIC_POLY(_ID)));

        polynomial_add(tmp_o, tmp_product, tmp_e, _ID + 1, tmp_summand, full_e_power2, &final_o,
                       C, (l * N_QUARTIC_POLY(_ID)), final_e);

        //Hint: Das hier funktioniert soweit gut für l1_Q2 (oft gedebuggt)
    }
}

void quartic_gf16v_madd2(uint8_t *C, const uint8_t *Av, unsigned A_pointer_index, char A_linear, const unsigned char *B,
                         unsigned B_pointer_index, unsigned B_offset, unsigned size_batch,
                         unsigned size_Bcolvec) {

    ///SHOULD BE DONE BETTER (WIP)--///
    int e_linear[2] = {2, 3}; // the structure of the sk-fields (do we need a constant factor?)
    unsigned char tmp_product[(N_QUARTIC_POLY(_ID) + 1) / 2]; // could be better calculated with i4.. in poly.c
    unsigned char tmp_summand[(_ID + 2) / 2]; //GF16, round up, one extra field for constant

    int tmp_e[15]; //size is too big..
    int final_e[15];

    int tmp_o = 0;
    int final_o = 0;

    //#define e_standard(n) (n+1) -> mal ersetzen
    int full_e_power2[15] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    int *e_A;
    int o_A;
    unsigned A_loop_offset;

    if (A_linear) {
        e_A = e_linear;
        o_A = 2;
        A_loop_offset = _ID;
    } else {
        e_A = full_e_power2;
        o_A = 6;
        A_loop_offset = N_QUARTIC_POLY(_ID);
    }
    ///--SHOULD BE DONE BETTER (WIP)///

    for (unsigned l = 0; l < size_batch * 2; l++) { // *2 for gf16 (size is in byte)
        //the inner loop of gf16vmadd

//        polynomial_print(o_A, &Av[(A_pointer_index) * size_batch * A_loop_offset], l * A_loop_offset, e_A, "l1_Q2:");
//        polynomial_print(2, &B[B_pointer_index * size_Bcolvec], B_offset * _ID, e_linear, "T1");

        polynomial_mul(o_A, &Av[(A_pointer_index) * size_batch * A_loop_offset], l * A_loop_offset, e_A, 2,
                       &B[B_pointer_index * size_Bcolvec],
                       B_offset * _ID,
                       e_linear, &tmp_o, tmp_product, 0, tmp_e);

//        polynomial_print(tmp_o, tmp_product, 0, tmp_e, "Produkt:");

        gf16_lin_poly_copy(tmp_summand, C, (l * N_QUARTIC_POLY(_ID)));

        polynomial_add(tmp_o, tmp_product, tmp_e, _ID + 1, tmp_summand, full_e_power2, &final_o,
                       C, (l * N_QUARTIC_POLY(_ID)), final_e);

//        polynomial_print(final_o,C,(l * N_QUARTIC_POLY(_ID)),final_e,"Written:");
    }
}

///
/// \param C
/// \param Av
/// \param A_pointer_index
/// \param A_linear
/// \param B
/// \param B_pointer_index
/// \param B_offset
/// \param size_batch
/// \param size_Bcolvec
/// @brief when C is linear over ID -> T4
void quartic_linear_gf16v_madd(uint8_t *C, const uint8_t *A, unsigned A_pointer_index, const unsigned char *B,
                               unsigned B_pointer_index, unsigned B_offset, unsigned size_batch,
                               unsigned size_Bcolvec) {

    ///SHOULD BE DONE BETTER (WIP)--///
    int e_linear[2] = {2, 3}; // the structure of the sk-fields (do we need a constant factor?)
    unsigned char tmp_product[(N_QUADRATIC_POLY(_ID) + 1) / 2]; // could be better calculated with i4.. in poly.c
    unsigned char tmp_summand[(_ID + 2) / 2]; //GF16, round up, one extra field for constant

    int tmp_e[15]; //size is too big..
    int final_e[15];

    int tmp_o = 0;
    int final_o = 0;

    ///--SHOULD BE DONE BETTER (WIP)///

    for (unsigned l = 0; l < size_batch * 2; l++) { // *2 for gf16 (size is in byte)
        //the inner loop of gf16vmadd
        polynomial_mul(2, &A[(A_pointer_index) * _ID * size_batch], l, e_linear, 2, &B[B_pointer_index * size_Bcolvec],
                       B_offset * _ID, // u.U l * _ID ????
                       e_linear, &tmp_o, tmp_product, 0, tmp_e);

        gf16_lin_poly_copy(tmp_summand, C, (l * _ID));

        polynomial_add(tmp_o, tmp_product, tmp_e, 2, tmp_summand, e_linear, &final_o,
                       C, (l * _ID), final_e);

        //Hint: Das hier funktioniert soweit gut für l1_Q2 (oft gedebuggt)
    }
}



