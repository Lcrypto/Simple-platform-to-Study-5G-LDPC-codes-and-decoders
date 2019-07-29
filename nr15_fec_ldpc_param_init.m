function [ldpc_param] = nr15_fec_ldpc_param_init(B,code_rate)

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

ldpc_param.set_iLS =    [2  4   8   16  32  64  128 256 ...
                         3  6   12  24  48  96  192 384 ...
                         5  10  20  40  80  160 320 ...
                         7  14  28  56  112 224 ...
                         9  18  36  72  144 288 ...
                         11 22  44  88  176 352 ...
                         13 26  52  104 208 ...
                         15 30  60  120 240];

ldpc_param.all_iLS =    [1  1   1   1   1   1   1   1 ...
                         2  2   2   2   2   2   2   2 ...
                         3  3   3   3   3   3   3 ...
                         4  4   4   4   4   4 ...
                         5  5   5   5   5   5 ...
                         6  6   6   6   6   6 ...
                         7  7   7   7   7  ...
                         8  8   8   8   8];

if B>3840
	A = B - 24;
else
	A = B - 16;
end
if (A<=292 || (A<=3824 && code_rate<=0.67) || code_rate<=0.25)
    base_graph_mode = 2;
    ldpc_param.K_cb = 3840;
else
    base_graph_mode = 1;
    ldpc_param.K_cb = 8448;
end
                     
set_iLS = ldpc_param.set_iLS;           
K_cb = ldpc_param.K_cb;

if (B <= K_cb)
    L = 0;
    C = 1;
    B1 = B;
else
    L = 24;
    C = ceil(B/(K_cb-L));
    B1 = B + C*L;
end

K_p = ceil(B1/C); % K_+ in the document
K_n = floor(B1/C);% K_- in the document

% % determine Z_c and K
if (base_graph_mode == 1)
    K_b = 22;
	Z_c = min(set_iLS(set_iLS*K_b - K_p >= 0));
    K = 22*Z_c;
else
    if (B > 640)
        K_b = 10;
    else
        if (B>560)
            K_b=9;
        else
            if (B>192)
                K_b=8;
            else
                K_b=6;
            end
        end
    end
    
	Z_c = min(set_iLS(set_iLS*K_b - K_p >= 0));
    
    if (isempty(Z_c))
        base_graph_mode = 1;
        K_b = 22;
        Z_c = min(set_iLS(set_iLS*K_b - K_p >= 0));
        K = 22*Z_c;
    else
        K = 10*Z_c;
    end
end

ldpc_param.K = K;
ldpc_param.B = B;
ldpc_param.K_b = K_b;
ldpc_param.K_p = K_p;
ldpc_param.K_n = K_n;
ldpc_param.C = C;
ldpc_param.L = L;
ldpc_param.iLS = ldpc_param.all_iLS(ldpc_param.set_iLS == Z_c);
ldpc_param.BG_sel = base_graph_mode;
ldpc_param.Z_c = Z_c;

if (base_graph_mode == 1)
    ldpc_param.N = 66*Z_c;
else
    ldpc_param.N = 50*Z_c;
end

[H_BG, H_row, H_col, H] = nr15_ldpc_matrix_init(ldpc_param);

ldpc_param.H_BG = H_BG;
ldpc_param.H_row = H_row;
ldpc_param.H_col = H_col;
ldpc_param.H = H;

% rate matching params
ldpc_param.code_rate = code_rate;
ldpc_param.k0 = 0;
ldpc_param.E = floor(K/code_rate);
ldpc_param.N_cb = ldpc_param.N;
