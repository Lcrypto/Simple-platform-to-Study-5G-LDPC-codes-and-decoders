
clear;
TBS = 100; %  the number of bits in the transport block (including CRC)
code_rate = 0.33;
snr = [0.3];
fer=0;
I_max = 50;
Frames=100;
%decode_mode = 'BP';
decode_mode = 'OMS';
beta = 0.6; % scale or offset parameter of the decoder
rng(100);
[ldpc_param] = ldpc_param_init(TBS,code_rate);



for i=1:Frames
tbs_bits = randi(2, 1, ldpc_param.B) - 1;

[cbs_bits] = nr15_ldpc_cbs(tbs_bits, ldpc_param);

[coded_bits, punctured_bits] = fec_ldpc_encoder(cbs_bits,ldpc_param);
[encoder_err] = fec_ldpc_check_enc(coded_bits, punctured_bits, ldpc_param.H); % just check that all ok



[rmed_bits] = fec_ldpc_rate_matching(coded_bits,ldpc_param);

mapped_sigs = 2.*rmed_bits-1;

received_sigs = awgn(mapped_sigs, snr+10*log10(code_rate/0.5));

sigma2 = 1/10^((snr+10*log10(code_rate/0.5))/10);
llr = -2*received_sigs./sigma2;

[dermed_llr] = fec_ldpc_rate_dematching(llr,ldpc_param);

[decoded_cbs_bits2] = fec_ldpc_decoder(dermed_llr, ldpc_param, I_max, decode_mode,beta);


output_bits = ldpc_decbs(decoded_cbs_bits2, ldpc_param);

err_bits = sum(abs(output_bits-tbs_bits));
ber = err_bits/TBS;
if (err_bits~=0)
fer=fer+1;

end
fer/i
end
% ldpc_param

% encoder_err

