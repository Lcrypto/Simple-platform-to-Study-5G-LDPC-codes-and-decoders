#include "ldpc.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void nr15_fec_ldpc_encoder_scb(int8_t* info_bits, nr15_ldpc_t *h, int8_t* coded_bits)
{
	int N = h->N;
	int K = h->K;
	int Z_c = h->Z_c;
	int16_t** H_BG = h->H_BG;

	int K_b, col_num;
	int Vij, Vp00, Vp10, Vp11, Vp20;

	int8_t *p0 = (int8_t*)malloc(sizeof(int8_t)*Z_c);
	int8_t *p1 = (int8_t*)malloc(sizeof(int8_t)*Z_c);
	int8_t *p2 = (int8_t*)malloc(sizeof(int8_t)*Z_c);
	int8_t *p3 = (int8_t*)malloc(sizeof(int8_t)*Z_c);
	int8_t *coded_temp = (int8_t*)malloc(sizeof(int8_t)*(N + 2 * Z_c));
	memset(p0, 0, sizeof(int8_t)*Z_c);
	memset(p1, 0, sizeof(int8_t)*Z_c);
	memset(p2, 0, sizeof(int8_t)*Z_c);
	memset(p3, 0, sizeof(int8_t)*Z_c);
	memset(coded_temp, 0, sizeof(int8_t)*(N + 2 * Z_c));
	memcpy(coded_temp, info_bits, K);

	if (h->BG_sel == 1)
	{
		K_b = 22;
		col_num = 68;
		Vp10 = H_BG[1][K_b];
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[0][j];
			if (Vij >= 0)
				circshift_xor(p1, coded_temp + j * Z_c, -Vij % Z_c, Z_c);
		}
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[1][j];
			if (Vij >= 0)
				circshift_xor(p2, coded_temp + j * Z_c, -Vij % Z_c, Z_c);
		}
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[2][j];
			if (Vij >= 0)
				circshift_xor(p3, coded_temp + j * Z_c, -Vij % Z_c, Z_c);
		}
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[3][j];
			if (Vij >= 0)
				circshift_xor(p0, coded_temp + j * Z_c, -(Vij - Vp10) % Z_c, Z_c);
		}

		circshift_xor(p0, p1, Vp10%Z_c, Z_c);
		circshift_xor(p0, p2, Vp10%Z_c, Z_c);
		circshift_xor(p0, p3, Vp10%Z_c, Z_c);
		memcpy(coded_temp + K_b * Z_c, p0, Z_c);

		Vp00 = H_BG[0][K_b];
		circshift_xor(p1, p0, -Vp00 % Z_c, Z_c);
		memcpy(coded_temp + (K_b + 1)*Z_c, p1, Z_c);

		Vp10 = H_BG[1][K_b];
		Vp11 = H_BG[1][K_b + 1];
		circshift_xor(p2, p0, -Vp10 % Z_c, Z_c);
		circshift_xor(p2, p1, -Vp11 % Z_c, Z_c);
		memcpy(coded_temp + (K_b + 2)*Z_c, p2, Z_c);

		circshift_xor(p3, p2, 0, Z_c);
		memcpy(coded_temp + (K_b + 3)*Z_c, p3, Z_c);

		for (int k = 4; k < col_num - K_b; k++)
		{
			for (int j = 0; j < K_b + k; j++)
			{
				Vij = H_BG[k][j];
				if (Vij > 0)
					circshift_xor(coded_temp + (K_b + k) * Z_c, coded_temp + j * Z_c, -Vij % Z_c, Z_c);
			}
		}
	}
	else
	{
		K_b = h->K_b;
		col_num = 52;

		Vp20 = H_BG[2][K_b];
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[0][j];
			if (Vij >= 0)
				circshift_xor(p1, coded_temp + Z_c * j, -Vij % Z_c, Z_c);
		}
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[1][j];
			if (Vij >= 0)
				circshift_xor(p2, coded_temp + Z_c * j, -Vij % Z_c, Z_c);
		}
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[2][j];
			if (Vij >= 0)
				circshift_xor(p3, coded_temp + Z_c * j, -Vij % Z_c, Z_c);
		}
		for (int j = 0; j < K_b; j++)
		{
			Vij = H_BG[3][j];
			if (Vij >= 0)
				circshift_xor(p0, coded_temp + Z_c * j, -(Vij - Vp20) % Z_c, Z_c);
		}

		circshift_xor(p0, p1, Vp20%Z_c, Z_c);
		circshift_xor(p0, p2, Vp20%Z_c, Z_c);
		circshift_xor(p0, p3, Vp20%Z_c, Z_c);
		memcpy(coded_temp + K_b * Z_c, p0, Z_c);

		Vp00 = H_BG[0][K_b];
		circshift_xor(p1, p0, -Vp00 % Z_c, Z_c);
		memcpy(coded_temp + (K_b + 1)*Z_c, p1, Z_c);

		circshift_xor(p2, p1, 0, Z_c);
		memcpy(coded_temp + (K_b + 2)*Z_c, p2, Z_c);

		Vp20 = H_BG[2][K_b];
		circshift_xor(p3, p0, -Vp20 % Z_c, Z_c);
		circshift_xor(p3, p2, 0, Z_c);
		memcpy(coded_temp + (K_b + 3)*Z_c, p3, Z_c);

		for (int k = 4; k < col_num - K_b; k++)
		{
			for (int j = 0; j < K_b + k; j++)
			{
				Vij = H_BG[k][j];
				if (Vij > 0)
					circshift_xor(coded_temp + (K_b + k) * Z_c, coded_temp + j * Z_c, -Vij % Z_c, Z_c);
			}
		}
	}

	memcpy(coded_bits, coded_temp + 2 * Z_c, N);

	free(coded_temp);
	free(p0);
	free(p1);
	free(p2);
	free(p3);
}

void circshift_xor(int8_t* p, int8_t* q, int num, int len)
{
	if (num < 0)
	{
		num = -num;
		for (int i = 0; i < num; i++)
			p[len + i - num] ^= q[i];
		for (int i = num; i < len; i++)
			p[i - num] ^= q[i];
	}
	else
	{
		for (int i = 0; i < num; i++)
			p[i] ^= q[len + i - num];
		for (int i = num; i < len; i++)
			p[i] ^= q[i - num];
	}
}