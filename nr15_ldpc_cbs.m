function [output_bits] = nr15_ldpc_cbs(input_bits, ldpc_param)

% Code block segmentation and code block CRC attachment 
% Ref: Section 5.2.1 in TS38.212

% Inputs and Ouputs:
% input_bits: binary input bits, size 1xB
% output_bits: binary output bits, size CxK

% Parameters and Variables:
% B: the length of TBS, 'B' in TS38.212
% K: the length of code blocks after CBS, 'K' in TS38.212
% code_rate: the code rate, 'R' in TS38.212
% set_iLS: Sets of LDPC lifting size, Table 5.3.2-1 in TS38.212
% K_cb: maximum code block size, 'K_cb' in TS38.212
% K_p, K_n: 'K_+'and 'K_-' in TS38.212
% K_b: 'K_b' in TS38.212
% L: CRC lenth, 'L' in TS38.212
% C: number of code block after CBS, 'C' in TS38.212
% B1: B plus the lenght of CRC, 'B'' in TS38.212

B = ldpc_param.B;
C = ldpc_param.C;
K = ldpc_param.K;
K_p = ldpc_param.K_p;
K_n = ldpc_param.K_n;
L = ldpc_param.L;

s = 0;
output_bits = -ones(C,K);

for r=0:C-1
    if r < mod(B,C)
        K1 = K_p;
    else
        K1 = K_n;
    end
    output_bits(r+1,1:K1-L) = input_bits(1,s+1:s+K1-L);
    if (C > 1)
        output_bits(r+1,K1-L+1:K1) = nr15_crc24b_encode(input_bits(1,s+1:s+K1-L));
    end
    s = s + K1-L;
end
