#include <mex.h>
#include <stdint.h>
#include <stdlib.h>

#include "ldpc.h"

/***************************************************************************
* Name:		mex_nr15_fec_ldpc_encoder_scb
*
* Function: ldpc scb encoder
*
* Input:	prhs[0]:	information bits
*			prhs[1]:	LDPC parameter array
*			prhs[2]:    base graph matrix
*
* Output:	plhs[0]:    coded bits
*			plhs[1]:    punctured bits
***************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* get prhs pointer */
    double *info_bits_D = mxGetPr(prhs[0]);
    double *ldpc_param_array_D = mxGetPr(prhs[1]);
    double *H_BG_D = mxGetPr(prhs[2]);

    /* format input */
    int32_t *ldpc_param_array =
        (int32_t *)malloc(sizeof(int32_t) * SIZE_OF_LDPC_ARRAY);
    for (int i = 0; i < SIZE_OF_LDPC_ARRAY; i++)
        ldpc_param_array[i] = ldpc_param_array_D[i];
    nr15_ldpc_t *ldpc_param =
        (nr15_ldpc_t *)malloc(sizeof(nr15_ldpc_t));
    trans_array_to_nr15_ldpc_t(ldpc_param, ldpc_param_array);
    ldpc_param->H_BG =
        (int16_t **)malloc(sizeof(int16_t *) * ldpc_param->row_hbg);
    ldpc_param->H_BG[0] =
        (int16_t *)malloc(sizeof(int16_t) * ldpc_param->row_hbg * ldpc_param->col_hbg);
    for (int i = 1; i < ldpc_param->row_hbg; i++)
        ldpc_param->H_BG[i] = ldpc_param->H_BG[i - 1] + ldpc_param->col_hbg;
    for (int i = 0; i < ldpc_param->row_hbg; i++)
        for (int j = 0; j < ldpc_param->col_hbg; j++)
            ldpc_param->H_BG[i][j] = H_BG_D[j * ldpc_param->row_hbg + i];
    int8_t *info_bits =
        (int8_t *)malloc(sizeof(int8_t) * ldpc_param->K);
    for (int i = 0; i < ldpc_param->K; i++)
        info_bits[i] = info_bits_D[i];

    /* encode */
    int8_t *coded_bits =
        (int8_t *)malloc(sizeof(int8_t) * ldpc_param->N);
    nr15_fec_ldpc_encoder_scb(info_bits, ldpc_param, coded_bits);

    /* get plhs pointer */
    plhs[0] = mxCreateDoubleMatrix(1, ldpc_param->N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, ldpc_param->Z_c * 2, mxREAL);

    /* format output */
    double *output = mxGetPr(plhs[0]);
    for (int i = 0; i < ldpc_param->N; i++)
        output[i] = coded_bits[i];
    double *punctured_bits = mxGetPr(plhs[1]);
    for (int i = 0; i < ldpc_param->Z_c * 2; i++)
        punctured_bits[i] = info_bits[i];

    free(ldpc_param_array);
    free_ldpc_param(ldpc_param);
    free(ldpc_param);
    free(info_bits);
    free(coded_bits);
}