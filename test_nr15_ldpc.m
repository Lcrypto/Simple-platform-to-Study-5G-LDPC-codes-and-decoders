% Test complete process of ldpc
clear;

B = 8448;
code_rate = 0.75;
snr = 3;
I_max = 50;
decode_mode = "BP";
beta = 0.6;

[ldpc_param] = nr15_fec_ldpc_param_init(B,code_rate);

tbs_bits = randi(2, 1, ldpc_param.B) - 1;

[cbs_bits] = nr15_ldpc_cbs(tbs_bits, ldpc_param);

[coded_bits, punctured_bits] = nr15_fec_ldpc_encoder_mex(cbs_bits,ldpc_param);

[encoder_err] = nr15_fec_ldpc_check_enc(coded_bits, punctured_bits, ldpc_param.H);

[rmed_bits] = nr15_fec_ldpc_rate_matching(coded_bits,ldpc_param);

mapped_sigs = 2.*rmed_bits-1;

received_sigs = awgn(mapped_sigs, snr+10*log10(code_rate/0.5));

sigma2 = 1/10^((snr+10*log10(code_rate/0.5))/10);
llr = -2*received_sigs./sigma2;

[dermed_llr] = nr15_fec_ldpc_rate_dematching(llr,ldpc_param);

[decoded_cbs_bits] = nr15_fec_ldpc_decoder_mex(dermed_llr, ldpc_param, I_max, decode_mode, beta);

output_bits = nr15_ldpc_decbs(decoded_cbs_bits, ldpc_param);

err_bits = sum(abs(output_bits-tbs_bits));
ber = err_bits/B

% ldpc_param

% encoder_err

