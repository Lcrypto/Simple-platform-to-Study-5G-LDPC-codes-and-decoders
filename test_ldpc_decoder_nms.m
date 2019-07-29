% Test nms decoder performance of different alpha
clear;

% simulation param
I_max = 100;
MAX_BLOCK = 1000;
MIN_ERR = 5000;
list_snr = -3:0.25:3;
list_indx = 1:length(list_snr);

list_coef = 0.64:-0.02:0.54;
list_indx_coef = 1:length(list_coef);

list_ber = zeros(length(list_coef),length(list_snr));
flag = ones(size(list_coef));

% init ldpc_param
B = 8448;
code_rate = 1/3;
fprintf("Start initializing ldpc_param...\n");
[ldpc_param] = nr15_fec_ldpc_param_init(B,code_rate);
ldpc_param_array = trans_ldpc_param_into_array(ldpc_param);
fprintf("Initializing ldpc_param done...\n");
K = ldpc_param.K;
Z_c = ldpc_param.Z_c;

err_bits = zeros(size(list_coef));
% start simulating
for indx_snr = list_indx
    snr = list_snr(indx_snr);
    
    sum_err = zeros(size(list_coef));
    sum_bits = 0;
    
    for indx_block = 1:MAX_BLOCK
        
        % get info_bits
        info_bits = randi(2, 1, K) - 1;
        
        % encode
        [coded_bits, punctured_bits] = nr15_fec_ldpc_encoder_scb(info_bits,ldpc_param);
        coded_bits_test = [punctured_bits, coded_bits];
        
        % map
        mapped_sigs = 2.*coded_bits_test-1;
        
        % pass channel
        received_sigs = awgn(mapped_sigs,snr);
        
        % get llr
        sigma2 = 1/10^(snr/10);
        llr = -2*received_sigs./sigma2;
        
        % decoder_ms
        for indx_coef = list_indx_coef
            if(flag(indx_coef))
                decoded_bits = mex_nr15_fec_ldpc_decoder_nms_punctured(...
                    llr(Z_c*2+1:end),ldpc_param_array,ldpc_param.H_col, ldpc_param.H_row,I_max,list_coef(indx_coef));
                err_bits(indx_coef) = sum(abs(decoded_bits - info_bits));
            else
                err_bits(indx_coef) = 0;
            end
            sum_err(indx_coef) = sum_err(indx_coef) + err_bits(indx_coef);
        end
        
        sum_bits = sum_bits + K;
        fprintf('No.%d\tSNR:%.2f(dB)\tBER:\t',indx_block,list_snr(indx_snr));
        for indx_coef = list_indx_coef
            fprintf('%.2f: %f\t',list_coef(indx_coef),err_bits(indx_coef)/K);
        end
        fprintf('\n');
        
        cnt_min_err = 0;
        for indx_coef = list_indx_coef
            if sum_err(indx_coef) > MIN_ERR
               cnt_min_err = cnt_min_err + 1; 
            end
        end
        if cnt_min_err >= length(list_indx_coef)
            break;
        end
    end
    
    % statistics
    for indx_coef = list_indx_coef
        list_ber(indx_coef,indx_snr) = sum_err(indx_coef)/sum_bits;
    end
    fprintf('=====SNR: %.2f(dB)\tBER:\t',list_snr(indx_snr));
    for indx_coef = list_indx_coef
        fprintf('%.2f:\t%f\t',list_coef(indx_coef),list_ber(indx_coef));
    end
    fprintf('=====\n');
    
    for indx_coef = list_indx_coef
        if sum_err(indx_coef) == 0
            flag(indx_coef) = 0;
        end
    end
    
    cnt_zero = 0;
    for indx_coef = list_indx_coef
        if sum_err(indx_coef) == 0
            cnt_zero = cnt_zero + 1;
        end
    end
    if cnt_zero >= length(list_indx)
        break;
    end
    
end

save('nms.mat','list_snr','list_ber');
