% Test performance of different decoders
clear;

% simulation param
I_max = 50;
MAX_BLOCK = 1000;
MIN_ERR = 10000;
list_snr = 0:0.25:3;
list_indx = 1:length(list_snr);
alpha = 0.797;
beta = 0.6;

% initial ber list and flag
list_ber_oms_pun = zeros(size(list_snr));
flag_oms_pun = 1;
list_ber_nms_pun = zeros(size(list_snr));
flag_nms_pun = 1;
list_ber_ms_pun = zeros(size(list_snr));
flag_ms_pun = 1;
list_ber_bp_pun = zeros(size(list_snr));
flag_bp_pun = 1;

% init ldpc_param
B = 8448;
code_rate = 0.75;
fprintf("Start initializing ldpc_param...\n");
[ldpc_param] = nr15_fec_ldpc_param_init(B,code_rate);
fprintf("Initializing ldpc_param done...\n");
K = ldpc_param.K;
Z_c = ldpc_param.Z_c;

% start simulating
for indx = list_indx
    snr = list_snr(indx);
    
    sum_err_oms_pun = 0;
    sum_err_nms_pun = 0;
    sum_err_ms_pun = 0;
    sum_err_bp_pun = 0;
    sum_bits = 0;
    
    for indx_block = 1:MAX_BLOCK
        
        % get tbs_bits
        tbs_bits = randi(2, 1, ldpc_param.B) - 1;
        
        % code block segmentation
        [cbs_bits] = nr15_ldpc_cbs(tbs_bits, ldpc_param);
        
        % encode
        [coded_bits, punctured_bits] = nr15_fec_ldpc_encoder_mex(cbs_bits,ldpc_param);
        
        % rate matching
        [rmed_bits] = nr15_fec_ldpc_rate_matching(coded_bits,ldpc_param);
        
        % map
        mapped_sigs = 2.*rmed_bits-1;
        
        % through channel
        received_sigs = awgn(mapped_sigs,snr+10*log10(code_rate/0.5));
        
        % get llr
        sigma2 = 1/10^((snr+10*log10(code_rate/0.5))/10);
        llr = -2*received_sigs./sigma2;
        
        % rate dematching
        dermed_llr = nr15_fec_ldpc_rate_dematching(llr,ldpc_param);
        
        % decoder_oms_pun
        if(flag_oms_pun)
            decoded_cbs_bits = nr15_fec_ldpc_decoder(dermed_llr,ldpc_param, I_max, "OMS", beta);
            decoded_bits = nr15_ldpc_decbs(decoded_cbs_bits, ldpc_param);
            err_bits_oms_pun = sum(abs(decoded_bits - tbs_bits));
        else
            err_bits_oms_pun = 0;
        end
        sum_err_oms_pun = sum_err_oms_pun + err_bits_oms_pun;
        
        % decoder_nms_pun
        if(flag_nms_pun)
            decoded_cbs_bits = nr15_fec_ldpc_decoder(dermed_llr,ldpc_param, I_max, "NMS", alpha);
            decoded_bits = nr15_ldpc_decbs(decoded_cbs_bits, ldpc_param);
            err_bits_nms_pun = sum(abs(decoded_bits - tbs_bits));
        else
            err_bits_nms_pun = 0;
        end
        sum_err_nms_pun = sum_err_nms_pun + err_bits_nms_pun;
        
        % decoder_ms_pun
        if(flag_ms_pun)
            decoded_cbs_bits = nr15_fec_ldpc_decoder(dermed_llr,ldpc_param, I_max, "MS", 0);
            decoded_bits = nr15_ldpc_decbs(decoded_cbs_bits, ldpc_param);
            err_bits_ms_pun = sum(abs(decoded_bits - tbs_bits));
        else
            err_bits_ms_pun = 0;
        end
        sum_err_ms_pun = sum_err_ms_pun + err_bits_ms_pun;
        
        % decoder_bp_pun
        if(flag_bp_pun)
            decoded_cbs_bits = nr15_fec_ldpc_decoder(dermed_llr,ldpc_param, I_max, "BP", 0);
            decoded_bits = nr15_ldpc_decbs(decoded_cbs_bits, ldpc_param);
            err_bits_bp_pun = sum(abs(decoded_bits - tbs_bits));
        else
            err_bits_bp_pun = 0;
        end
        sum_err_bp_pun = sum_err_bp_pun + err_bits_bp_pun;
        
        % statistics of one block
        sum_bits = sum_bits + K;
        fprintf('No.%d\tSNR:%.2f(dB)\tBER:\toms: %f\tnms: %f\tms: %f\tbp: %f\n',...
            indx_block,...
            list_snr(indx),...
            err_bits_oms_pun/K,...
            err_bits_nms_pun/K,...
            err_bits_ms_pun/K,...
            err_bits_bp_pun/K);
        if (sum_err_oms_pun > MIN_ERR || flag_oms_pun == 0) &&...
           (sum_err_nms_pun > MIN_ERR || flag_nms_pun == 0) &&...
           (sum_err_ms_pun > MIN_ERR || flag_ms_pun == 0) &&...
           (sum_err_bp_pun > MIN_ERR || flag_bp_pun == 0)
            break;
        end
    end
    
    % statistics
    list_ber_oms_pun(indx) = sum_err_oms_pun/sum_bits;
    list_ber_nms_pun(indx) = sum_err_nms_pun/sum_bits;
    list_ber_ms_pun(indx) = sum_err_ms_pun/sum_bits;
    list_ber_bp_pun(indx) = sum_err_bp_pun/sum_bits;
    fprintf('=====SNR: %.2f(dB)\tBER:\toms:%f\tnms:%f\tms:%f\tbp:%f=====\n',...
        list_snr(indx),...
        list_ber_oms_pun(indx),...
        list_ber_nms_pun(indx),...
        list_ber_ms_pun(indx),...
        list_ber_bp_pun(indx));
    if sum_err_oms_pun == 0
        flag_oms_pun = 0;
    end
    if sum_err_nms_pun == 0
        flag_nms_pun = 0;
    end
    if sum_err_ms_pun == 0
        flag_ms_pun = 0;
    end
    if sum_err_bp_pun == 0
        flag_bp_pun = 0;
    end
    if sum_err_oms_pun == 0 && sum_err_nms_pun == 0 && sum_err_ms_pun == 0 && sum_err_bp_pun == 0
        break;
    end
end
