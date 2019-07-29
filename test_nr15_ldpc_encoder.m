% Test ldpc encoder
clear;

B = 8448;
code_rate = 0.5;

[ldpc_param] = nr15_fec_ldpc_param_init(B,code_rate);

tbs_bits = randi(2, 1, ldpc_param.B) - 1;

[cbs_bits] = nr15_ldpc_cbs(tbs_bits, ldpc_param);

[coded_bits, punctured_bits] = nr15_fec_ldpc_encoder(cbs_bits,ldpc_param);
[coded_bits_mex, punctured_bits_mex] = nr15_fec_ldpc_encoder_mex(cbs_bits,ldpc_param);

[encoder_err] = nr15_fec_ldpc_check_enc(coded_bits, punctured_bits, ldpc_param.H);
[encoder_err_mex] = nr15_fec_ldpc_check_enc(coded_bits_mex, punctured_bits_mex, ldpc_param.H);

ldpc_param

encoder_err

