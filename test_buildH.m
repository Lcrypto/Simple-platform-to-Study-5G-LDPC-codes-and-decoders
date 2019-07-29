clear;

[ldpc_param] = nr15_ldpc_init();
ldpc_param.Z_c = 224;
ldpc_param.base_graph_mode = 1;

Z_c = ldpc_param.Z_c;
base_graph_mode = ldpc_param.base_graph_mode;
i_LS = ldpc_param.all_iLS(ldpc_param.set_iLS == Z_c);

if (ldpc_param.base_graph_mode == 1)
    row = 46;
    col = 68;
else
    row = 42;
    col = 52;
end

IZc = eye(Z_c);
H_BG = load('BG1_iLS4');
H = zeros(row*ldpc_param.Z_c ,col*ldpc_param.Z_c);

for i = 1:row
	for j = 1:col
        Vij = H_BG(i,j);
        if Vij >= 0
            H((i-1)*Z_c+1:(i-1)*Z_c+Z_c,(j-1)*Z_c+1:(j-1)*Z_c+Z_c) = circshift(IZc,[0 mod(Vij,Z_c)]);
        end
	end
end

% [max(sum(H,2)) max(sum(H,1))]

N = size(H, 2);
K = N - size(H, 1);

A = zeros(N,max(sum(H,1))+1);          % Row i contains the indexes of the check nodes in which the variable i is involved
B = zeros(N-K,max(sum(H,2))+1);        % Row i contains the indexes of the variable nodes involved in the check i

% for vn=1:N
%     tmph_col = find(H(:,vn));
%     tmp_col_weight = length(tmph_col);
%     A(vn,1) = tmp_col_weight;
%     A(vn,2:tmp_col_weight+1) = tmph_col;
% end
% 
% for cn=1:(N-K)
%     tmph_row = find(H(cn,:));
%     tmp_row_weight = length(tmph_row);
%     B(cn,1) = tmp_row_weight;
%     B(cn,2:tmp_row_weight+1) = tmph_row;
% end

B = 4928; 
ldpc_param.code_rate = 0.5;
input_bits = randi(2, 1, B) - 1;

[output_bits, ldpc_param.Z_c, ldpc_param.base_graph_mode] = nr15_ldpc_cbs(input_bits, ldpc_param);

K_b = 22;

coded_bits = reshape(zeros(1, N), Z_c, []);

coded_bits(:,1:K_b) = reshape(output_bits, Z_c, []);

p0 = zeros(Z_c,1);
for j = 1:K_b
    Vij = H_BG(1,j);
    if (Vij >=0)
        p0 = p0 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
	Vij = H_BG(2,j);
    if (Vij >=0)
        p0 = p0 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
	Vij = H_BG(3,j);
    if (Vij >=0)
        p0 = p0 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
	Vij = H_BG(4,j);
    if (Vij >=0)
        p0 = p0 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
end

Vp10 = H_BG(2,K_b+1);
shift_p0_10 = p0; 
p0 = circshift(shift_p0_10, mod(Vp10,Z_c));
p0 = mod(p0,2);
Vp00 = H_BG(1,K_b+1);
shift_p0_00 = circshift(p0,-Vp00);
coded_bits(:,K_b+1) = p0;

p1 = shift_p0_00;
for j = 1:K_b
    Vij = H_BG(1,j);
    if (Vij >=0)
        p1 = p1 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
end
p1 = mod(p1,2);
coded_bits(:,K_b+2) = p1;

p2 = shift_p0_00;
for j = 1:K_b
    Vij = H_BG(3,j);
    if (Vij >=0)
        p2 = p2 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
	Vij = H_BG(4,j);
    if (Vij >=0)
        p2 = p2 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
end
p2 = mod(p2,2);
coded_bits(:,K_b+3) = p2;

p3 = shift_p0_00;
for j = 1:K_b
    Vij = H_BG(4,j);
    if (Vij >=0)
        p3 = p3 + circshift(coded_bits(:,j),-mod(Vij,Z_c));
    end
end
p3 = mod(p3,2);
coded_bits(:,K_b+4) = p3;

for k=5:col-K_b
    tmp_pk = zeros(Z_c,1);
    for j = 1:K_b+k-1
        Vij = H_BG(k,j);
        if (Vij >=0)
            tmp_pk = tmp_pk + circshift(coded_bits(:,j),-mod(Vij,Z_c));
        end
    end
    tmp_pk = mod(tmp_pk,2);
    coded_bits(:,K_b+k) = tmp_pk;
end

coded_bits_vec = reshape(coded_bits, [], 1);

