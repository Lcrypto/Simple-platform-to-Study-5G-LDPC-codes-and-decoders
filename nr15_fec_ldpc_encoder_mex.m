% ldpc encoder with mex function
function [coded_bits, punctured_bits] = nr15_fec_ldpc_encoder_mex(cbs_bits,ldpc_param)

N = ldpc_param.N;
Z_c = ldpc_param.Z_c;
C = ldpc_param.C;

coded_bits = zeros(C, N);
punctured_bits = zeros(C, 2*Z_c);
ldpc_param_array = trans_ldpc_param_into_array(ldpc_param);

for k=1:C
    [coded_bits(k,:), punctured_bits(k,:)] = ...
        mex_nr15_fec_ldpc_encoder_scb(cbs_bits(k,:),ldpc_param_array, ldpc_param.H_BG);
end

end