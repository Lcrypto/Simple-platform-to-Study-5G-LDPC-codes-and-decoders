#ifndef LDPC_H
#define LDPC_H

#include <stdint.h>
#include <stdlib.h>

#define ILS_NUM 51
#define MAX_SET_ILS 384
#define SIZE_OF_LDPC_ARRAY 17

typedef struct variable_node
{
	int8_t degree;	// number of connective check nodes
	int16_t *index;	// index of connective check nodes
	float* message;	// message from vn to cn
	int8_t pointer;	// pointer to the message that will be used
	float judge_message;
} variable_node;

typedef struct check_node
{
	int8_t degree;	// number of connective variable nodes
	int16_t *index;	// index of connective variable nodes
	float* message;	// message from cn to vn
	int8_t pointer;	// pointer to the message that will be used
} check_node;

typedef struct nr15_ldpc_t
{
	int B;					// the length of TBS, 'B' in TS38.212
	int code_rate;			// code rate*1024
	int K_b;				// block number of message bits, 'K_b' in TS38.212
	int K_p;				// 'K_-' in TS38.212
	int K_n;				// 'K_+' in TS38.212
	int C;					// number of code block after CBS, 'C' in TS38.212
	int L;					// CRC lenth, 'L' in TS38.212
	int iLS;				// LDPC lifting size
	int BG_sel;				// number of base graph
	int Z_c;				// block size, 'Z_c' in TS38.212
	int K;					// message bits length, 'K' in TS38.212
	int N;					// matrix length, 'N' in TS38.212
	int col_hbg;			// column number of base graph
	int row_hbg;			// row number of base graph
	int max_col_weight;		// max column weight of H
	int max_row_weight;		// max row weight of H
	int K_cb;				// maximum code block size, 'K_cb' in TS38.212
	int16_t **H_BG;			// base graph
	variable_node* vn_list;	// variable node list
	check_node* cn_list;	// check node list
} nr15_ldpc_t;

/*************************************************************************************/
/*                            Declare LDPC initial functions                         */
/*************************************************************************************/

/* initial LDPC parameter */
void nr15_fec_ldpc_param_init(nr15_ldpc_t *h, int B, double code_rate);

/* initial base graph matrix */
void nr15_ldpc_matrix_init(nr15_ldpc_t *h);

/* initial variable node list and check node list */
void nr15_ldpc_node_list_init(nr15_ldpc_t *h);

/* transform ldpc_param structure into parameter array */
void trans_nr15_ldpc_t_to_array(nr15_ldpc_t *h, int32_t* array);

/* transform parameter array into ldpc_param structure*/
void trans_array_to_nr15_ldpc_t(nr15_ldpc_t *h, int32_t* array);

void free_ldpc_param_with_node(nr15_ldpc_t *h);

void free_ldpc_param(nr15_ldpc_t *h);

void free_variable_node(variable_node* vn);

void free_check_node(check_node* cn);

/*************************************************************************************/
/*                            Declare LDPC encode functions                          */
/*************************************************************************************/

/*
* Name:		nr15_fec_ldpc_encoder_scb
*
* Input:	info_bits:	information array
*			h:			LDPC type structure
*
* Output:	coded_bits:	coded array
*/
void nr15_fec_ldpc_encoder_scb(int8_t* info_bits, nr15_ldpc_t *h, int8_t* coded_bits);

void circshift_xor(int8_t* p, int8_t* q, int num, int len);

/*************************************************************************************/
/*                            Declare LDPC decode functions                          */
/*************************************************************************************/

/*
* Name:		nr15_fec_ldpc_decoder_oms_punctured
*
* Input:	llr:			received llr message
*			ldpc_arg:		LDPC type structure
*			I_max:			max iteration
*			beta:			offset parameter
*
* Output:	decoded_bits:	decoded array
*/
void nr15_fec_ldpc_decoder_oms_punctured(float *llr, nr15_ldpc_t *ldpc_arg, int I_max, float beta, int8_t *decoded_bits);

/*
* Name:		nr15_fec_ldpc_decoder_nms_punctured
*
* Input:	llr:			received llr message
*			ldpc_arg:		LDPC type structure
*			I_max:			max iteration
*			alpha:			normalize parameter
*
* Output:	decoded_bits:	decoded array
*/
void nr15_fec_ldpc_decoder_nms_punctured(float *llr, nr15_ldpc_t *ldpc_arg, int I_max, float alpha, int8_t *decoded_bits);

/*
* Name:		nr15_fec_ldpc_decoder_ms_punctured
*
* Input:	llr:			received llr message
*			ldpc_arg:		LDPC type structure
*			I_max:			max iteration
*
* Output:	decoded_bits:	decoded array
*/
void nr15_fec_ldpc_decoder_ms_punctured(float *llr, nr15_ldpc_t *ldpc_arg, int I_max, int8_t *decoded_bits);

/*
* Name:		nr15_fec_ldpc_decoder_bp_punctured
*
* Input:	llr:			received llr message
*			ldpc_arg:		LDPC type structure
*			I_max:			max iteration
*
* Output:	decoded_bits:	decoded array
*/
void nr15_fec_ldpc_decoder_bp_punctured(float *llr, nr15_ldpc_t *ldpc_arg, int I_max, int8_t *decoded_bits);

int check_decoded_bits(int8_t* whole_decoded_bits, nr15_ldpc_t* h);

#endif