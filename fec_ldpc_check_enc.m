function [encoder_err] = fec_ldpc_check_enc(coded_bits, punctured_bits, H)

tmp_bits = coded_bits;
tmp_bits(tmp_bits<0) = 0;

encoder_err = sum(abs(mod(H*[punctured_bits tmp_bits].',2)));

