% simplified rate dematching
function [dermed_llr] = rate_dematching(llr,ldpc_param)

E = ldpc_param.E;
% K = ldpc_param.K;
N = ldpc_param.N;
C = ldpc_param.C;
dermed_llr = zeros(C,N);

for c = 1:C
    if E<=N
        dermed_llr(c,1:E) = llr(c,1:E);
    else
        dermed_llr(c,1:N) = llr(c,1:N);
        dermed_llr(c,1:E-N) = dermed_llr(c,1:E-N) + llr(c,E+1:end);
    end
end

end