function [output_bits] = nr15_ldpc_decbs(input_bits, ldpc_param)

C = ldpc_param.C;
K = ldpc_param.K;
B = ldpc_param.B;
K_p = ldpc_param.K_p;
K_n = ldpc_param.K_n;
L = ldpc_param.L;

rem = mod(B, C);

output_bits = zeros(1,B);

num = 1;
for c = 1:C
    if c <= rem
        output_bits(num:num+K_p-L-1) = input_bits(c,1:K_p-L);
        num = num + K_p - L;
    else
        output_bits(num:num+K_n-L-1) = input_bits(c,1:K_n-L);
        num = num + K_n - L;
    end
end

end