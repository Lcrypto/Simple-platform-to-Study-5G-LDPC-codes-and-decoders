% bp decoder without 2Z_c puncturing
function [decoded_bits,check_bits] = nr15_fec_ldpc_decoder_bp(llr,ldpc_param,I_max)

N = size(ldpc_param.H,2);
M = size(ldpc_param.H,1);
K = ldpc_param.K;
H_row = ldpc_param.H_row;
H_col = ldpc_param.H_col;
H = ldpc_param.H;

vn = zeros(M,N);
cn = zeros(M,N);
w = zeros(N,0);

% Initialize variable nodes
for n = 1:N
    for m = H_col(n,2:(H_col(n,1)+1))
        vn(m,n) = llr(n);
    end
end

for i = 1:I_max
    % Update the check nodes
    for m = 1:M
        for n = H_row(m,2:(H_row(m,1)+1))
            array_n1 = H_row(m,2:(H_row(m,1)+1));
            array_n1(array_n1==n) = [];
%             cn(m,n) = min(abs(vn(m,array_n1)))*prod(sign(vn(m,array_n1)));
            cn(m,n)=2*atanh(prod(tanh(vn(m,array_n1)./2)));
        end
    end
    
    % Update the variable nodes
    for n = 1:N
        for m = H_col(n,2:(H_col(n,1)+1))
            array_m1 = H_col(n,2:(H_col(n,1)+1));
            array_m1(array_m1==m) = [];
            vn(m,n) = llr(n) + sum(cn(array_m1,n));
        end
    end
    
    % Apply a hard decision
    for n =1:N
        if llr(n)+sum(cn(H_col(n,2:(H_col(n,1)+1)),n))>=0
            w(n) = 0;
        else
            w(n) = 1;
        end
    end
    
    % Check decision
    if sum(mod(w*H',2))==0
        break;
    end
    
end
decoded_bits = w(1:K);
check_bits = w(K+1:end);

end