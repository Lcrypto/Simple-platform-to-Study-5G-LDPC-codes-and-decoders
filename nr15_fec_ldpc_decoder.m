% ldpc decoder with m function
function decoded_bits= nr15_fec_ldpc_decoder(llr,ldpc_param, I_max, decode_mode, coef)

K = ldpc_param.K;
C = ldpc_param.C;

decoded_bits = zeros(C, K);

% choose decode mode
if(decode_mode == "NMS")
    for k=1:C
        decoded_bits(k,:) = nr15_fec_ldpc_decoder_nms_punctured(llr(k,:), ldpc_param, I_max,coef);
    end
elseif(decode_mode == "OMS")
    for k=1:C
        decoded_bits(k,:) = nr15_fec_ldpc_decoder_oms_punctured(llr(k,:), ldpc_param, I_max,coef);
    end
elseif(decode_mode == "MS")
    for k=1:C
        decoded_bits(k,:) = nr15_fec_ldpc_decoder_ms_punctured(llr(k,:), ldpc_param, I_max);
    end
else
    for k=1:C
        decoded_bits(k,:) = nr15_fec_ldpc_decoder_bp_punctured(llr(k,:), ldpc_param, I_max);
    end
end
end