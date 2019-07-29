function [H] = nr15_ldpc_H_init(ldpc_param)

BG_sel = ldpc_param.BG_sel;
iLS = ldpc_param.iLS;
Z_c = ldpc_param.Z_c;

H_BG = load(['BG',num2str(BG_sel),'_','iLS',num2str(iLS)]);

[H_BG_row_num, H_BG_col_num] = size(H_BG);

IZc = eye(Z_c);
H = zeros(H_BG_row_num*ldpc_param.Z_c ,H_BG_col_num*ldpc_param.Z_c);

for i = 1:H_BG_row_num
	for j = 1:H_BG_col_num
        Vij = H_BG(i,j);
        if Vij >= 0
            H((i-1)*Z_c+1:(i-1)*Z_c+Z_c,(j-1)*Z_c+1:(j-1)*Z_c+Z_c) = circshift(IZc,[0 mod(Vij,Z_c)]);
        end
	end
end
