% rate match by nr15 method
function [rmed_bits] = nr15_fec_ldpc_rate_matching(coded_bits,ldpc_param)

E = ldpc_param.E;
k0 = ldpc_param.k0;
N = ldpc_param.N;
N_cb = ldpc_param.N_cb;
C = ldpc_param.C;
rmed_bits = zeros(C,E);

for c = 1:C
    k = 1;
    j = 1;
    while(k<=E)
        if mod(k0+j-1,N_cb)+1<=N
            rmed_bits(c,k) = coded_bits(c,mod(k0+j-1,N_cb)+1);
            k = k + 1;
        end
        j = j + 1;
    end
end

end