% rate dematch by nr15 method
function [dermed_llr] = nr15_fec_ldpc_rate_dematching(llr,ldpc_param)

E = ldpc_param.E;
k0 = ldpc_param.k0;
N = ldpc_param.N;
N_cb = ldpc_param.N_cb;
C = ldpc_param.C;
dermed_llr = zeros(C,N);

for c = 1:C
    k = 1;
    j = 1;
    while(k<=E)
        if mod(k0+j-1,N_cb)+1<=N
            dermed_llr(c,mod(k0+j-1,N_cb)+1) = dermed_llr(c,mod(k0+j-1,N_cb)+1) + llr(c,k);
            k = k + 1;
        end
        j = j + 1;
    end
end

end