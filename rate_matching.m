% simplified rate matching
function [rmed_bits] = rate_matching(coded_bits,ldpc_param)

E = ldpc_param.E;
N = ldpc_param.N;
C = ldpc_param.C;
rmed_bits = zeros(C,E);

for c = 1:C
    if E<=N
        rmed_bits(c,1:E) = coded_bits(c,1:E);
    else
        rmed_bits(c,1:N) = coded_bits(c,1:N);
        rmed_bits(c,N+1:end) = coded_bits(c,1:E-N);
    end
end

end