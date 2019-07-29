#include "ldpc.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int set_iLS[ILS_NUM] = { 2, 4, 8, 16, 32, 64, 128, 256,
							   3, 6, 12, 24, 48, 96, 192, 384,
							   5, 10, 20, 40, 80, 160, 320,
							   7, 14, 28, 56, 112, 224,
							   9, 18, 36, 72, 144, 288,
							   11, 22, 44, 88, 176, 352,
							   13, 26, 52, 104, 208,
							   15, 30, 60, 120, 240 };

static int all_iLS[ILS_NUM] = { 1, 1, 1, 1, 1, 1, 1, 1,
							   2, 2, 2, 2, 2, 2, 2, 2,
							   3, 3, 3, 3, 3, 3, 3,
							   4, 4, 4, 4, 4, 4,
							   5, 5, 5, 5, 5, 5,
							   6, 6, 6, 6, 6, 6,
							   7, 7, 7, 7, 7,
							   8, 8, 8, 8, 8 };

void nr15_fec_ldpc_param_init(nr15_ldpc_t *h, int B, double code_rate)
{
	int A, B1, L, K_p, K_n, K_b, Z_c, K_cb, C, K, N;
	int base_graph_mode, iLS, row_hbg, col_hbg;

	if (B > 3840)
		A = B - 24;
	else
		A = B - 16;

	if (A <= 292 ||
		(A <= 3824 && code_rate <= 0.67) ||
		code_rate <= 0.25)
	{
		base_graph_mode = 2;
		h->K_cb = 3840;
	}
	else
	{
		base_graph_mode = 1;
		h->K_cb = 8448;
	}

	K_cb = h->K_cb;

	if (B <= K_cb)
	{
		L = 0;
		C = 1;
		B1 = B;
	}
	else
	{
		L = 24;
		C = (B - 1) / (K_cb - L) + 1;
		B1 = B + C * L;
	}

	K_p = (B1 - 1) / C + 1;
	K_n = B1 / C;

	if (base_graph_mode == 1)
	{
		K_b = 22;
		Z_c = MAX_SET_ILS + 1;
		for (int i = 0; i < ILS_NUM; i++)
			if (set_iLS[i] * K_b - K_p >= 0 && set_iLS[i] < Z_c)
			{
				Z_c = set_iLS[i];
				iLS = all_iLS[i];
			}

		K = 22 * Z_c;
	}
	else
	{
		if (B > 640)
			K_b = 10;
		else if (B > 560)
			K_b = 9;
		else if (B > 192)
			K_b = 8;
		else
			K_b = 6;
		Z_c = MAX_SET_ILS + 1;
		for (int i = 0; i < ILS_NUM; i++)
			if (set_iLS[i] * K_b - K_p >= 0 && set_iLS[i] < Z_c)
			{
				Z_c = set_iLS[i];
				iLS = all_iLS[i];
			}

		if (Z_c == MAX_SET_ILS + 1)
		{
			base_graph_mode = 1;
			K_b = 22;
			Z_c = MAX_SET_ILS + 1;
			for (int i = 0; i < ILS_NUM; i++)
				if (set_iLS[i] * K_b - K_p >= 0 && set_iLS[i] < Z_c)
				{
					Z_c = set_iLS[i];
					iLS = all_iLS[i];
				}
			K = 22 * Z_c;
		}
		else
			K = 10 * Z_c;
	}

	if (base_graph_mode == 1)
	{
		row_hbg = 46;
		col_hbg = 68;
		N = 66 * Z_c;
	}
	else
	{
		row_hbg = 42;
		col_hbg = 52;
		N = 50 * Z_c;
	}

	h->K = K;
	h->B = B;
	h->K_b = K_b;
	h->K_p = K_p;
	h->K_n = K_n;
	h->C = C;
	h->L = L;
	h->iLS = iLS;
	h->BG_sel = base_graph_mode;
	h->Z_c = Z_c;
	h->N = N;
	h->row_hbg = row_hbg;
	h->col_hbg = col_hbg;
	h->code_rate = code_rate * 1024;

	// initial H_BG
	nr15_ldpc_matrix_init(h);
	nr15_ldpc_node_list_init(h);

}

void nr15_ldpc_matrix_init(nr15_ldpc_t *h)
{
	int row_hbg = h->row_hbg;
	int col_hbg = h->col_hbg;
	int BG_sel = h->BG_sel;
	int iLS = h->iLS;

	h->H_BG = (int16_t **)malloc(sizeof(int16_t*) * row_hbg);
	h->H_BG[0] = (int16_t *)malloc(sizeof(int16_t) * row_hbg * col_hbg);
	for (int i = 1; i < row_hbg; i++)
		h->H_BG[i] = h->H_BG[i - 1] + col_hbg;

	int16_t** H_BG = h->H_BG;

	char hbg_name[9] = "BG0_iLS0";
	hbg_name[2] = hbg_name[2] + (char)BG_sel;
	hbg_name[7] = hbg_name[7] + (char)iLS;
	FILE *fbg;
	if ((fbg = fopen(hbg_name, "r")) == NULL)
	{
		printf("Can not open %s !\n", hbg_name);
		getch();
		exit(0);
	}
	for (int i = 0; i < row_hbg; i++)
		for (int j = 0; j < col_hbg; j++)
			fscanf(fbg, "%d", &H_BG[i][j]);
	fclose(fbg);
}

void nr15_ldpc_node_list_init(nr15_ldpc_t *h)
{
	int i, j, m, n;
	int N = h->N + 2 * h->Z_c;
	int M = N - h->K;
	h->cn_list = (check_node*)malloc(sizeof(check_node)*M);
	for (m = 0; m < M; m++)
		h->cn_list[m].degree = 0;
	h->vn_list = (variable_node*)malloc(sizeof(variable_node)*N);
	for (n = 0; n < N; n++)
		h->vn_list[n].degree = 0;

	for (i = 0; i < h->row_hbg; i++)
		for (j = 0; j < h->col_hbg; j++)
		{
			if (h->H_BG[i][j] >= 0)
			{
				for (m = i * h->Z_c; m < h->Z_c*(i + 1); m++)
				{
					n = j * h->Z_c + (m + h->H_BG[i][j]) % h->Z_c;
					//h->cn_list[m].index[h->cn_list[m].degree] = n;
					h->cn_list[m].degree++;
					//h->vn_list[n].index[h->vn_list[n].degree] = m;
					h->vn_list[n].degree++;
				}
			}
		}

	for (n = 0; n < N; n++)
	{
		h->vn_list[n].index = (int16_t*)malloc(sizeof(int16_t)*h->vn_list[n].degree);
		h->vn_list[n].message = (float*)malloc(sizeof(float)*h->vn_list[n].degree);
		//if (n < 100)
		//	printf("vn[%d]:%d\n", n, h->vn_list[n].degree);
		h->vn_list[n].degree = 0;
	}

	for (m = 0; m < M; m++)
	{
		h->cn_list[m].index = (int16_t*)malloc(sizeof(int16_t)*h->cn_list[m].degree);
		h->cn_list[m].message = (float*)malloc(sizeof(float)*h->cn_list[m].degree);
		//if (m < 100)
		//	printf("cn[%d]:%d\n", m, h->cn_list[m].degree);
		h->cn_list[m].degree = 0;
	}

	for (i = 0; i < h->row_hbg; i++)
		for (j = 0; j < h->col_hbg; j++)
		{
			if (h->H_BG[i][j] >= 0)
			{
				for (m = i * h->Z_c; m < h->Z_c*(i + 1); m++)
				{
					n = j * h->Z_c + (m + h->H_BG[i][j]) % h->Z_c;
					h->cn_list[m].index[h->cn_list[m].degree] = n;
					h->cn_list[m].degree++;
					h->vn_list[n].index[h->vn_list[n].degree] = m;
					h->vn_list[n].degree++;
				}
			}
		}
	h->max_col_weight = h->vn_list[0].degree;
	h->max_row_weight = h->cn_list[0].degree;
	//for (int m = 0; m < 100; m++)
	//{
	//	printf("cn[%d]:", m);
	//	for (int i = 0; i < h->cn_list[m].degree; i++)
	//		printf("%d ", h->cn_list[m].index[i]);
	//	printf("\n");
	//}
	//for (int n = 0; n < 100; n++)
	//{
	//	printf("vn[%d]:", n);
	//	for (int i = 0; i < h->vn_list[n].degree; i++)
	//		printf("%d ", h->vn_list[n].index[i]);
	//	printf("\n");
	//}

}

/*
fec_ldpc_param[0]: B;
fec_ldpc_param[1]: code_rate; code_rate = fec_ldpc_param[1]/100;
fec_ldpc_param[2]: K_b;
fec_ldpc_param[3]: K_p;
fec_ldpc_param[4]: K_n;
fec_ldpc_param[5]: C;
fec_ldpc_param[6]: L;
fec_ldpc_param[7]: iLS;
fec_ldpc_param[8]: BG_sel;
fec_ldpc_param[9]: Z_c;
fec_ldpc_param[10]: K;
fec_ldpc_param[11]: N;
fec_ldpc_param[12]: h_col_num
fec_ldpc_param[13]: h_row_num
fec_ldpc_param[14]: max_col_weight
fec_ldpc_param[15]: max_row_weight
fec_ldpc_param[16]: K_cb
*/
void trans_nr15_ldpc_t_to_array(nr15_ldpc_t *h, int32_t* fec_ldpc_param)
{
	fec_ldpc_param[0] = h->B;
	fec_ldpc_param[1] = h->code_rate;
	fec_ldpc_param[2] = h->K_b;
	fec_ldpc_param[3] = h->K_p;
	fec_ldpc_param[4] = h->K_n;
	fec_ldpc_param[5] = h->C;
	fec_ldpc_param[6] = h->L;
	fec_ldpc_param[7] = h->iLS;
	fec_ldpc_param[8] = h->BG_sel;
	fec_ldpc_param[9] = h->Z_c;
	fec_ldpc_param[10] = h->K;
	fec_ldpc_param[11] = h->N;
	fec_ldpc_param[12] = h->col_hbg;
	fec_ldpc_param[13] = h->row_hbg;
	fec_ldpc_param[14] = h->max_col_weight;
	fec_ldpc_param[15] = h->max_row_weight;
	fec_ldpc_param[16] = h->K_cb;
}

void trans_array_to_nr15_ldpc_t(nr15_ldpc_t *h, int32_t* fec_ldpc_param)
{
	h->B = fec_ldpc_param[0];
	h->code_rate = fec_ldpc_param[1];
	h->K_b = fec_ldpc_param[2];
	h->K_p = fec_ldpc_param[3];
	h->K_n = fec_ldpc_param[4];
	h->C = fec_ldpc_param[5];
	h->L = fec_ldpc_param[6];
	h->iLS = fec_ldpc_param[7];
	h->BG_sel = fec_ldpc_param[8];
	h->Z_c = fec_ldpc_param[9];
	h->K = fec_ldpc_param[10];
	h->N = fec_ldpc_param[11];
	h->col_hbg = fec_ldpc_param[12];
	h->row_hbg = fec_ldpc_param[13];
	h->max_col_weight = fec_ldpc_param[14];
	h->max_row_weight = fec_ldpc_param[15];
	h->K_cb = fec_ldpc_param[16];
}

void free_ldpc_param(nr15_ldpc_t *h)
{
	free(h->H_BG[0]);
	free(h->H_BG);
}

void free_ldpc_param_with_node(nr15_ldpc_t *h)
{
	int N = h->N + 2 * h->Z_c;
	int M = N - h->K;
	for (int n = 0; n < N; n++)
		free_variable_node(&(h->vn_list[n]));
	free(h->vn_list);
	for (int m = 0; m < M; m++)
		free_check_node(&(h->cn_list[m]));
	free(h->cn_list);
}

void free_variable_node(variable_node* vn)
{
	free(vn->index);
	free(vn->message);
}

void free_check_node(check_node* cn)
{
	free(cn->index);
	free(cn->message);
}