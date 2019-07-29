function [H_BG, H_row, H_col, H] = nr15_ldpc_matrix_init(ldpc_param)

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

[H_row_num, H_col_num] = size(H);

H_col = zeros(H_col_num,max(sum(H,1))+1);          % Row i contains the indexes of the check nodes in which the variable i is involved
H_row = zeros(H_row_num,max(sum(H,2))+1);        % Row i contains the indexes of the variable nodes involved in the check i

for j=1:H_col_num
    tmph_col = find(H(:,j));
    tmp_col_weight = length(tmph_col);
    H_col(j,1) = tmp_col_weight;
    H_col(j,2:tmp_col_weight+1) = tmph_col;
end

for i=1:H_row_num
    tmph_row = find(H(i,:));
    tmp_row_weight = length(tmph_row);
    H_row(i,1) = tmp_row_weight;
    H_row(i,2:tmp_row_weight+1) = tmph_row;
end
