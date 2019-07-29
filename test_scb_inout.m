clear;

B = 3872;
code_rate = 0.5;

[ldpc_param] = nr15_fec_ldpc_param_init(B,code_rate);

cbs_info_bits=load('enc_scb_in.txt');

[coded_bits, punctured_bits] = nr15_fec_ldpc_encoder_scb(cbs_info_bits.',ldpc_param);

matlab_cbs_coded_bits=[punctured_bits, coded_bits].';

cbs_coded_bits=load('enc_scb_out.txt');

err=(abs(cbs_coded_bits - matlab_cbs_coded_bits));

sum(err)