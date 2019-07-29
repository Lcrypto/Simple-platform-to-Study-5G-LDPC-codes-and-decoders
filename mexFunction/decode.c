#include "ldpc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_LLR 1E20

double excepted_prod(double* array, int indx, int len);

void nr15_fec_ldpc_decoder_oms_punctured(float *llr, nr15_ldpc_t *h, int I_max, float beta, int8_t *decoded_bits)
{
	int i, j, m, n, pointer;
	int temp_sign, indx_min;
	float temp_llr, min_llr, submin_llr, sum_llr;
	int N = h->N + 2 * h->Z_c;
	int M = N - h->K;

	/* initial nodes' pointers */
	check_node *cn_list = h->cn_list;
	variable_node *vn_list = h->vn_list;
	for (n = 0; n < N; n++)
		vn_list[n].pointer = 0;
	for (m = 0; m < M; m++)
		cn_list[m].pointer = 0;

	/* initial temporary memory */
	float *whole_llr = (float*)malloc(sizeof(float)*N);
	memset(whole_llr, 0, 2 * h->Z_c * sizeof(float));
	memcpy(whole_llr + 2 * h->Z_c, llr, h->N * sizeof(float));
	int8_t *whole_decoded_bits = (int8_t*)malloc(sizeof(int8_t)*N);
	int8_t *sign_list =	(int8_t*)malloc(sizeof(int8_t)*h->col_hbg);

	/* initialize vn_list */
	for (int n = 0; n < N; n++)
	{
		for (i = 0; i < vn_list[n].degree; i++)
		{
			vn_list[n].message[i] = whole_llr[n];
		}
	}

	int indx_i;
	for (indx_i = 0; indx_i < I_max; indx_i++)
	{
		/* upadate the check nodes */
		for (m = 0; m < M; m++)
		{
			temp_sign = 1;
			min_llr = MAX_LLR;
			submin_llr = MAX_LLR;
			for (i = 0; i < cn_list[m].degree; i++)
			{
				n = cn_list[m].index[i];
				pointer = vn_list[n].pointer;
				temp_llr = vn_list[n].message[pointer];
				vn_list[n].pointer = (pointer + 1) % vn_list[n].degree;
				/* calculate sign */
				if (temp_llr < 0)
				{
					temp_sign = -temp_sign;
					temp_llr = -temp_llr;
					sign_list[i] = -1;
				}
				else
					sign_list[i] = 1;
				/* calculate min LLR */
				if (temp_llr < min_llr)
				{
					submin_llr = min_llr;
					min_llr = temp_llr;
					indx_min = i;
				}
				else if (temp_llr < submin_llr)
					submin_llr = temp_llr;
			}
			min_llr -= beta;
			submin_llr -= beta;
			min_llr = min_llr > 0 ? min_llr : 0;
			submin_llr = submin_llr > 0 ? submin_llr : 0;
			for (i = 0; i < cn_list[m].degree; i++)
			{
				if (i == indx_min)
					cn_list[m].message[i] = temp_sign * sign_list[i] * submin_llr;
				else
					cn_list[m].message[i] = temp_sign * sign_list[i] * min_llr;
			}
		}

		/* update the variable nodes */
		for (n = 0; n < N; n++)
		{
			sum_llr = 0.0;
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				sum_llr += cn_list[m].message[pointer];
			}
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				vn_list[n].message[j] = whole_llr[n] + sum_llr - cn_list[m].message[pointer];
				cn_list[m].pointer = (pointer + 1) % cn_list[m].degree;
			}

			/* hard judgement */
			whole_decoded_bits[n] = (sum_llr + whole_llr[n]) > 0 ? 0 : 1;
		}

		if (check_decoded_bits(whole_decoded_bits, h))
			break;
	}
	memcpy(decoded_bits, whole_decoded_bits, h->K);
	free(whole_decoded_bits);
	free(whole_llr);
	free(sign_list);
}

void nr15_fec_ldpc_decoder_nms_punctured(float *llr, nr15_ldpc_t *h, int I_max, float alpha, int8_t *decoded_bits)
{
	int i, j, m, n, pointer;
	int temp_sign, indx_min;
	float temp_llr, min_llr, submin_llr, sum_llr;
	int N = h->N + 2 * h->Z_c;
	int M = N - h->K;

	/* initial nodes' pointers */
	check_node *cn_list = h->cn_list;
	variable_node *vn_list = h->vn_list;
	for (n = 0; n < N; n++)
		vn_list[n].pointer = 0;
	for (m = 0; m < M; m++)
		cn_list[m].pointer = 0;

	/* initial temporary memory */
	float *whole_llr = (float*)malloc(sizeof(float)*N);
	memset(whole_llr, 0, 2 * h->Z_c * sizeof(float));
	memcpy(whole_llr + 2 * h->Z_c, llr, h->N * sizeof(float));
	int8_t *whole_decoded_bits = (int8_t*)malloc(sizeof(int8_t)*N);
	int8_t *sign_list = (int8_t*)malloc(sizeof(int8_t)*h->col_hbg);

	/* Initialize vn_list */
	for (int n = 0; n < N; n++)
	{
		for (i = 0; i < vn_list[n].degree; i++)
		{
			vn_list[n].message[i] = whole_llr[n];
		}
	}

	int indx_i;
	for (indx_i = 0; indx_i < I_max; indx_i++)
	{
		/* upadate the check nodes */
		for (m = 0; m < M; m++)
		{
			temp_sign = 1;
			min_llr = MAX_LLR;
			submin_llr = MAX_LLR;
			for (i = 0; i < cn_list[m].degree; i++)
			{
				n = cn_list[m].index[i];
				pointer = vn_list[n].pointer;
				temp_llr = vn_list[n].message[pointer];
				vn_list[n].pointer = (pointer + 1) % vn_list[n].degree;
				/* calculate sign */
				if (temp_llr < 0)
				{
					temp_sign = -temp_sign;
					temp_llr = -temp_llr;
					sign_list[i] = -1;
				}
				else
					sign_list[i] = 1;
				/* calculate min LLR */
				if (temp_llr < min_llr)
				{
					submin_llr = min_llr;
					min_llr = temp_llr;
					indx_min = i;
				}
				else if (temp_llr < submin_llr)
					submin_llr = temp_llr;
			}
			for (i = 0; i < cn_list[m].degree; i++)
			{
				if (i == indx_min)
					cn_list[m].message[i] = alpha * temp_sign * sign_list[i] * submin_llr;
				else
					cn_list[m].message[i] = alpha * temp_sign * sign_list[i] * min_llr;
			}
		}

		/* update the variable nodes */
		for (n = 0; n < N; n++)
		{
			sum_llr = 0.0;
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				sum_llr += cn_list[m].message[pointer];
			}
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				vn_list[n].message[j] = whole_llr[n] + sum_llr - cn_list[m].message[pointer];
				cn_list[m].pointer = (pointer + 1) % cn_list[m].degree;
			}

			/* hard judgement */
			whole_decoded_bits[n] = (sum_llr + whole_llr[n]) > 0 ? 0 : 1;
		}

		if (check_decoded_bits(whole_decoded_bits, h))
			break;
	}
	memcpy(decoded_bits, whole_decoded_bits, h->K);
	free(whole_decoded_bits);
	free(whole_llr);
	free(sign_list);
}

void nr15_fec_ldpc_decoder_ms_punctured(float *llr, nr15_ldpc_t *h, int I_max, int8_t *decoded_bits)
{
	int i, j, m, n, pointer;
	int temp_sign, indx_min;
	float temp_llr, min_llr, submin_llr, sum_llr;
	int N = h->N + 2 * h->Z_c;
	int M = N - h->K;

	/* initial nodes' pointers */
	check_node *cn_list = h->cn_list;
	variable_node *vn_list = h->vn_list;
	for (n = 0; n < N; n++)
		vn_list[n].pointer = 0;
	for (m = 0; m < M; m++)
		cn_list[m].pointer = 0;

	/* initial temporary memory */
	float *whole_llr = (float*)malloc(sizeof(float)*N);
	memset(whole_llr, 0, 2 * h->Z_c * sizeof(float));
	memcpy(whole_llr + 2 * h->Z_c, llr, h->N * sizeof(float));
	int8_t *whole_decoded_bits = (int8_t*)malloc(sizeof(int8_t)*N);
	int8_t *sign_list = (int8_t*)malloc(sizeof(int8_t)*h->col_hbg);

	/* Initialize vn_list */
	for (int n = 0; n < N; n++)
	{
		for (i = 0; i < vn_list[n].degree; i++)
		{
			vn_list[n].message[i] = whole_llr[n];
		}
	}

	int indx_i;
	for (indx_i = 0; indx_i < I_max; indx_i++)
	{
		/* upadate the check nodes */
		for (m = 0; m < M; m++)
		{
			temp_sign = 1;
			min_llr = MAX_LLR;
			submin_llr = MAX_LLR;
			for (i = 0; i < cn_list[m].degree; i++)
			{
				n = cn_list[m].index[i];
				pointer = vn_list[n].pointer;
				temp_llr = vn_list[n].message[pointer];
				vn_list[n].pointer = (pointer + 1) % vn_list[n].degree;
				/* calculate sign */
				if (temp_llr < 0)
				{
					temp_sign = -temp_sign;
					temp_llr = -temp_llr;
					sign_list[i] = -1;
				}
				else
					sign_list[i] = 1;
				/* calculate min LLR */
				if (temp_llr < min_llr)
				{
					submin_llr = min_llr;
					min_llr = temp_llr;
					indx_min = i;
				}
				else if (temp_llr < submin_llr)
					submin_llr = temp_llr;
			}
			for (i = 0; i < cn_list[m].degree; i++)
			{
				if (i == indx_min)
					cn_list[m].message[i] = temp_sign * sign_list[i] * submin_llr;
				else
					cn_list[m].message[i] = temp_sign * sign_list[i] * min_llr;
			}
		}

		/* update the variable nodes */
		for (n = 0; n < N; n++)
		{
			sum_llr = 0.0;
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				sum_llr += cn_list[m].message[pointer];
			}
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				vn_list[n].message[j] = whole_llr[n] + sum_llr - cn_list[m].message[pointer];
				cn_list[m].pointer = (pointer + 1) % cn_list[m].degree;
			}

			/* hard judgement */
			whole_decoded_bits[n] = (sum_llr + whole_llr[n]) > 0 ? 0 : 1;
		}

		if (check_decoded_bits(whole_decoded_bits, h))
			break;
	}
	memcpy(decoded_bits, whole_decoded_bits, h->K);
	free(whole_decoded_bits);
	free(whole_llr);
	free(sign_list);
}

void nr15_fec_ldpc_decoder_bp_punctured(float *llr, nr15_ldpc_t *h, int I_max, int8_t *decoded_bits)
{
	int i, j, m, n, pointer, temp_sign;
	double temp_tanh;
	float temp_llr, sum_llr;
	int N = h->N + 2 * h->Z_c;
	int M = N - h->K;

	/* initial nodes' pointers */
	check_node *cn_list = h->cn_list;
	variable_node *vn_list = h->vn_list;
	for (n = 0; n < N; n++)
		vn_list[n].pointer = 0;
	for (m = 0; m < M; m++)
		cn_list[m].pointer = 0;

	/* initial temporary memory */
	float *whole_llr = (float*)malloc(sizeof(float)*N);
	memset(whole_llr, 0, 2 * h->Z_c * sizeof(float));
	memcpy(whole_llr + 2 * h->Z_c, llr, h->N * sizeof(float));
	int8_t *whole_decoded_bits = (int8_t*)malloc(sizeof(int8_t)*N);
	double *tanh_list = (double*)malloc(sizeof(double)*h->col_hbg);
	int8_t *sign_list = (int8_t*)malloc(sizeof(int8_t)*h->col_hbg);

	/* Initialize vn_list */
	for (int n = 0; n < N; n++)
	{
		for (i = 0; i < vn_list[n].degree; i++)
		{
			vn_list[n].message[i] = whole_llr[n];
		}
	}

	int indx_i;
	for (indx_i = 0; indx_i < I_max; indx_i++)
	{
		/* upadate the check nodes */
		for (m = 0; m < M; m++)
		{
			temp_sign = 1;
			for (i = 0; i < cn_list[m].degree; i++)
			{
				n = cn_list[m].index[i];
				pointer = vn_list[n].pointer;
				temp_llr = vn_list[n].message[pointer];
				vn_list[n].pointer = (pointer + 1) % vn_list[n].degree;
				/* calculate sign */
				if (temp_llr < 0)
				{
					temp_sign = -temp_sign;
					temp_llr = -temp_llr;
					sign_list[i] = -1;
				}
				else
					sign_list[i] = 1;
				/* calculate tanh */
				tanh_list[i] = tanh(temp_llr / 2);
			}
			for (i = 0; i < cn_list[m].degree; i++)
			{
				double temp = excepted_prod(tanh_list, i, cn_list[m].degree);
				if (temp >= 1)
					cn_list[m].message[i] = MAX_LLR * temp_sign*sign_list[i];
				else
					cn_list[m].message[i] = 2 * atanh(excepted_prod(tanh_list, i, cn_list[m].degree))*temp_sign*sign_list[i];
			}
		}

		/* update the variable nodes */
		for (n = 0; n < N; n++)
		{
			sum_llr = 0.0;
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				sum_llr += cn_list[m].message[pointer];
			}
			for (j = 0; j < vn_list[n].degree; j++)
			{
				m = vn_list[n].index[j];
				pointer = cn_list[m].pointer;
				vn_list[n].message[j] = whole_llr[n] + sum_llr - cn_list[m].message[pointer];
				cn_list[m].pointer = (pointer + 1) % cn_list[m].degree;
			}

			/* hard judgement */
			whole_decoded_bits[n] = (sum_llr + whole_llr[n]) > 0 ? 0 : 1;
		}

		if (check_decoded_bits(whole_decoded_bits, h))
			break;
	}
	memcpy(decoded_bits, whole_decoded_bits, h->K);
	free(whole_decoded_bits);
	free(whole_llr);
	free(tanh_list);
	free(sign_list);
}

int check_decoded_bits(int8_t* whole_decoded_bits, nr15_ldpc_t* h)
{
	int m, n, i;
	int N = h->N + 2 * h->Z_c;
	int M = N - h->K;
	check_node* cn_list = h->cn_list;

	int8_t temp;
	for (m = 0; m < M; m++)
	{
		temp = 0;
		for (i = 0; i < cn_list[m].degree; i++)
		{
			n = cn_list[m].index[i];
			temp ^= whole_decoded_bits[n];
		}
		if (temp != 0)
			return 0;
	}
	return 1;
}

double excepted_prod(double* array, int indx, int len)
{
	double ret = 1.0;
	for (int i = 0; i < len; i++)
		if (i != indx)
			ret *= array[i];
	return ret;
}