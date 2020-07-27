///  @file parallel_matrix_op.c
///  @brief the standard implementations for functions in parallel_matrix_op.h
///
///  the standard implementations for functions in parallel_matrix_op.h
///


#include "blas_comm.h"
#include "blas.h"
#include "rainbow_blas.h"

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


void quartic_UpperTrianglize_gf16(unsigned char *btriC, const unsigned char *bA, unsigned Awidth, unsigned size_batch) {

    unsigned char *runningC = btriC;
    unsigned Aheight = Awidth;

    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < i; j++) {
            unsigned idx = idx_of_trimat(j, i, Aheight);
            for (unsigned k = 0; k < size_batch * 2; k++) { //*2 because GF16

                polynomial_add(btriC + idx * size_batch * N_QUARTIC_POLY, N_QUARTIC_POLY * k, 3,
                               bA + size_batch * (i * Awidth + j) * N_QUARTIC_POLY, N_QUARTIC_POLY * k,
                               N_CUBIC_POLY, _full_e_power2);
            }
            //non quartic: gf256v_add( btriC + idx*size_batch , bA + size_batch*(i*Awidth+j) , size_batch );
        }
        for (unsigned l = 0; l < size_batch * (Aheight - i) * 2; l++) {
            polynomial_add(runningC, N_QUARTIC_POLY * l, 3, bA + size_batch * (i * Awidth + i) * N_QUARTIC_POLY,
                           l * N_QUARTIC_POLY, N_CUBIC_POLY, _full_e_power2);
        }
        //non quartic: gf256v_add(runningC, bA + size_batch * (i * Awidth + i), size_batch * (Aheight - i)); /// ATTENTION GF256
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
    //TODO: this function is not correct!
    unsigned Awidth = Bheight;
    unsigned Aheight = Awidth;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (k < i) continue;
                quartic_gf16v_madd(bC, 2, 4, &btriA[(k - i) * size_batch * N_LINEAR_POLY], 1, 1,
                                   &B[j * size_Bcolvec * _ID], k, 0, 0, size_batch);
//                gf16v_madd(bC, &btriA[(k - i) * size_batch], gf16v_get_ele(&B[j * size_Bcolvec], k), size_batch);
            }
            bC += (size_batch * N_QUARTIC_POLY);
        }
        btriA += (Aheight - i) * N_LINEAR_POLY * size_batch;
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

void quartic_batch_trimatTr_madd_gf16(unsigned char *bC, const unsigned char *btriA, const unsigned char *B,
                                      unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch) {
    unsigned Aheight = Bheight;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (i < k) continue;
                quartic_gf16v_madd(bC, 2, 4, &btriA[size_batch * (idx_of_trimat(k, i, Aheight)) * N_LINEAR_POLY], 1, 1,
                                   &B[j * size_Bcolvec * _ID], k, 0, 0,
                                   size_batch);
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
                              unsigned Awidth, const unsigned char *bB, unsigned Bwidth, unsigned size_batch) {
    unsigned Atr_height = Awidth;
    unsigned Atr_width = Aheight;
    for (unsigned i = 0; i < Atr_height; i++) {
        for (unsigned j = 0; j < Atr_width; j++) {
            quartic_gf16v_madd(bC, 3, 4, &bB[j * Bwidth * size_batch * N_QUARTIC_POLY], 2, 4,
                               &A_to_tr[size_Acolvec * i * _ID], j, 0, 0, size_batch * Bwidth);
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
                quartic_gf16v_madd(bC, 2, 4, &bA[size_batch * (i + k * Aheight) * N_LINEAR_POLY], 1, 1,
                                   &B[j * size_Bcolvec * _ID], k, 0, 0,
                                   size_batch);
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
                quartic_gf16v_madd(bC, 2, 4, &bA[k * size_batch * N_LINEAR_POLY], 1, 1, &B[j * size_Bcolvec * _ID], k,
                                   0,
                                   0, size_batch);
                //gf16v_madd( bC , & bA[ k*size_batch ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch * N_QUARTIC_POLY;
        }
        bA += (Awidth) * size_batch * N_LINEAR_POLY;
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

void quartic_gf16v_madd(uint8_t *C, unsigned C_grade, unsigned C_structure_grade, const uint8_t *A, unsigned A_grade,
                        unsigned A_strucutre_grade, const unsigned char *B, unsigned B_offset, unsigned B_grade,
                        unsigned B_structure_grade, unsigned size_batch) {
    //todo: NOT CORRECT FOR ODD id-size
    unsigned A_loop_offset = _grade_n_poly_terms(A_strucutre_grade);
    unsigned B_loop_offset = _grade_n_poly_terms(B_structure_grade);
    unsigned C_loop_offset = _grade_n_poly_terms(C_structure_grade);

    ///temporal variables to konow the structure of the product:
    unsigned char tmp_product[(N_QUARTIC_POLY + 5) / 2];
    unsigned tmp_e[N_QUARTIC_POLY + N_QUADRATIC_POLY]; //don't ask
    unsigned tmp_o = 0;

    for (unsigned l = 0; l < size_batch * 2; l++) { // *2 for gf16 (size is in byte)

        polynomial_mul(A, l * A_loop_offset, A_grade, B, B_offset * B_loop_offset, B_grade, tmp_product, &tmp_o, tmp_e);

        polynomial_add(C, l * C_loop_offset, C_grade, tmp_product, 0, tmp_o, tmp_e);
    }
}

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

    unsigned s_t_size = sizeof(((usk_t *) 0)->s1) + sizeof(((usk_t *) 0)->t1) + sizeof(((usk_t *) 0)->t4) +
                        sizeof(((usk_t *) 0)->t3);

    for (unsigned i = 0; i < s_t_size * 2; i++) {
        gf16v_set_ele(usk, i, polynomial_value(_ID, msk, i * _ID, _lin_e_power2, id));
    }
    usk += s_t_size;
    msk += s_t_size * _ID;

    for (unsigned i = 0; i < (sizeof(usk_t) - LEN_SKSEED - (s_t_size)) * 2; i++) {//+t4_size+t3_size)) * 2; i++) {
        gf16v_set_ele(usk, i, polynomial_value(N_LINEAR_POLY, msk, i * N_LINEAR_POLY, _full_e_power2, id));
    }
}


