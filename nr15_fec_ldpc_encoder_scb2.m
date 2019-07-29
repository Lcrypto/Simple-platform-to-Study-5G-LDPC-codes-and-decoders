function [coded_bits, punctured_bits] = nr15_fec_ldpc_encoder_scb(info_bits,ldpc_param)

N = ldpc_param.N;
K = ldpc_param.K;
Z_c = ldpc_param.Z_c;
H_BG = ldpc_param.H_BG;

coded_bits = zeros(1,N);

info_bits_matrix = reshape(info_bits, Z_c, []);
info_bits_matrix(info_bits_matrix<0) = 0;
punctured_bits = reshape(info_bits_matrix(:,1:2),1,[]);

if ldpc_param.BG_sel == 1
    K_b = 22; 
    col_num = 68;
    
    coded_bits_matrix = zeros(Z_c, col_num);
    coded_bits_matrix(:, 1:K/Z_c) = info_bits_matrix;

    % parity_bits_matrix = zeros(Z_c, (N-K)/Z_c+2);
    coded_bits(1, 1:K-2*Z_c) = info_bits(1,2*Z_c+1:K);
    
    p0 = zeros(Z_c,1);
    for j = 1:K_b
        Vij = H_BG(1,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
        Vij = H_BG(2,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
        Vij = H_BG(3,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
        Vij = H_BG(4,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end

    Vp10 = H_BG(2,K_b+1);
    shift_p0_10 = p0; 
    p0 = circshift(shift_p0_10, mod(Vp10,Z_c));
    p0 = mod(p0,2);
    Vp00 = H_BG(1,K_b+1);
    shift_p0_00 = circshift(p0,-Vp00);
    coded_bits_matrix(:,K_b+1) = p0;

    p1 = shift_p0_00;
    for j = 1:K_b
        Vij = H_BG(1,j);
        if (Vij >=0)
            p1 = p1 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end
    p1 = mod(p1,2);
    coded_bits_matrix(:,K_b+2) = p1;

    p2 = shift_p0_00;
    for j = 1:K_b
        Vij = H_BG(3,j);
        if (Vij >=0)
            p2 = p2 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
        Vij = H_BG(4,j);
        if (Vij >=0)
            p2 = p2 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end
    p2 = mod(p2,2);
    coded_bits_matrix(:,K_b+3) = p2;

    p3 = shift_p0_00;
    for j = 1:K_b
        Vij = H_BG(4,j);
        if (Vij >=0)
            p3 = p3 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end
    p3 = mod(p3,2);
    coded_bits_matrix(:,K_b+4) = p3;

    for k=5:col_num-K_b
        tmp_pk = zeros(Z_c,1);
        for j = 1:K_b+k-1
            Vij = H_BG(k,j);
            if (Vij >=0)
                tmp_pk = tmp_pk + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
            end
        end
        tmp_pk = mod(tmp_pk,2);
        coded_bits_matrix(:,K_b+k) = tmp_pk;
    end
else
    K_b = 10; % Here K_b is not the value defined in TS38.212
    col_num = 52;
    
    coded_bits_matrix = zeros(Z_c, col_num);
    coded_bits_matrix(:, 1:K/Z_c) = info_bits_matrix;
    
    coded_bits(1, 1:K-2*Z_c) = info_bits(1,2*Z_c+1:K);
    
	p0 = zeros(Z_c,1);
    for j = 1:K_b
        Vij = H_BG(1,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
        Vij = H_BG(2,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
        Vij = H_BG(3,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
        Vij = H_BG(4,j);
        if (Vij >=0)
            p0 = p0 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end

    Vp20 = H_BG(3,K_b+1);
    shift_p0_20 = p0; 
    p0 = circshift(shift_p0_20, mod(Vp20,Z_c));
    p0 = mod(p0,2);
    Vp00 = H_BG(1,K_b+1);
    shift_p0_00 = circshift(p0,-Vp00);
    Vp30 = H_BG(4,K_b+1);
    shift_p0_30 = circshift(p0,-Vp30);
    coded_bits_matrix(:,K_b+1) = p0;

    p1 = shift_p0_00;
    for j = 1:K_b
        Vij = H_BG(1,j);
        if (Vij >=0)
            p1 = p1 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end
    p1 = mod(p1,2);
    coded_bits_matrix(:,K_b+2) = p1;
    
    p2 = p1;
    for j = 1:K_b
        Vij = H_BG(2,j);
        if (Vij >=0)
            p2 = p2 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end
    p2 = mod(p2,2);
    coded_bits_matrix(:,K_b+3) = p2;

    p3 = shift_p0_30;
    for j = 1:K_b
        Vij = H_BG(4,j);
        if (Vij >=0)
            p3 = p3 + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
        end
    end
    p3 = mod(p3,2);
    coded_bits_matrix(:,K_b+4) = p3;

    for k=5:col_num-K_b
        tmp_pk = zeros(Z_c,1);
        for j = 1:K_b+k-1
            Vij = H_BG(k,j);
            if (Vij >=0)
                tmp_pk = tmp_pk + circshift(coded_bits_matrix(:,j),-mod(Vij,Z_c));
            end
        end
        tmp_pk = mod(tmp_pk,2);
        coded_bits_matrix(:,K_b+k) = tmp_pk;
    end
end

% coded_bits_vec = reshape(coded_bits, [], 1);
coded_bits(1, K-2*Z_c+1:N) = reshape(coded_bits_matrix(:,K_b+1:end),1,[]);

