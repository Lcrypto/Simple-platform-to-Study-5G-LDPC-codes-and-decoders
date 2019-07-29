% oms decoder with cutting H matrix
function [decoded_bits,check_bits] = nr15_fec_ldpc_decoder_oms_punctured(llr,ldpc_param,I_max,beta)

K = ldpc_param.K;
% N = size(ldpc_param.H,2);
N = ldpc_param.E + 2*ldpc_param.Z_c;
M = N-K;
H_row = ldpc_param.H_row;
H_col = ldpc_param.H_col;
H = ldpc_param.H;
Z_c = ldpc_param.Z_c;

vn = zeros(M,N);
cn = zeros(M,N);
vi = [zeros(1,2*Z_c) llr];
w = zeros(N,0);

% Initialize variable nodes
for n = 1:N
    array_m = H_col(n,2:(H_col(n,1)+1));
    array_m(array_m>M) = [];
    for m = array_m
        vn(m,n) = vi(n);
    end
end

for i = 1:I_max
    % Update the check nodes
    for m = 1:M
        array_n = H_row(m,2:(H_row(m,1)+1));
        array_n(array_n>N) = [];
        for n = array_n
            array_n1 = H_row(m,2:(H_row(m,1)+1));
            array_n1(array_n1==n) = [];
            array_n1(array_n1>N) = [];
            cn(m,n) = max(min(abs(vn(m,array_n1)))-beta,0)*prod(sign(vn(m,array_n1)));
        end
    end
    
    % Update the variable nodes
    for n = 1:N
        array_m = H_col(n,2:(H_col(n,1)+1));
        array_m(array_m>M) = [];
        for m = array_m
            array_m1 = H_col(n,2:(H_col(n,1)+1));
            array_m1(array_m1==m) = [];
            array_m1(array_m1>M) = [];
            vn(m,n) = vi(n) + sum(cn(array_m1,n));
        end
    end
    
    % Apply a hard decision
    for n =1:N
        array_m = H_col(n,2:(H_col(n,1)+1));
        array_m(array_m>M) = [];
        if vi(n)+sum(cn(array_m,n))>=0
            w(n) = 0;
        else
            w(n) = 1;
        end
    end
    
    % Check decision
    if sum(mod(w*H(1:M,1:N)',2))==0
        break;
    end
    
end
decoded_bits = w(1:K);
check_bits = w(K+1:end);

end