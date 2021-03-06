% ldpc encoder with m function
function [coded_bits, punctured_bits] = fec_ldpc_encoder(cbs_bits,ldpc_param)

N = ldpc_param.N;
Z_c = ldpc_param.Z_c;
C = ldpc_param.C;

coded_bits = zeros(C, N);
punctured_bits = zeros(C, 2*Z_c);

for k=1:C
    [coded_bits(k,:), punctured_bits(k,:)] = fec_ldpc_encoder_scb(cbs_bits(k,:),ldpc_param);
%     [coded_bits2, punctured_bits2] = fec_ldpc_encoder_scb2(cbs_bits(k,:),ldpc_param);
end

end