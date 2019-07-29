load('ber_r34_K3872.mat');

figure;

% semilogy([1.00, 1.25],[1.721e-3, 1.560e-7],'kd-');
% semilogy([1.00, 1.25, 1.50],[1.084e-1, 5.050e-2, 1.038e-6],'kd-');
% semilogy([1.00 1.25 1.50 1.75 2.00 2.25 2.50],...
%     [1.219e-1 1.132e-1 1.022e-1 8.779e-2 6.246e-2 9.630e-3 8.816e-6],'kd-');
% semilogy([1.00 1.25 1.50 1.75 2.00 2.25 2.50 2.75],...
%     [1.212e-01, 1.120e-01, 9.964e-02, 8.288e-02, 5.231e-02, 7.933e-03, 7.800e-05, 1.086e-08],'kd-');
semilogy([1.00 1.25 1.50 1.75 2.00 2.25 2.50 2.75 3.00],...
    [1.220e-01, 1.131e-01, 1.017e-01, 8.686e-02, 5.963e-02, 1.654e-02, 4.249e-04, 1.209e-06, 3.482e-10],'kd-');
hold on;
semilogy(list_snr,list_ber_bp_pun,'bs-');
semilogy(list_snr,list_ber_ms_pun,'go-');
semilogy(list_snr,list_ber_nms_pun,'mx-');
semilogy(list_snr,list_ber_oms_pun,'r+-');

legend('reference',...
    'BP',...
    'MS',...
    'NMS,\alpha=0.797',...
    'OMS,\beta=0.6'...
    );

axis([0 3 10e-5 1]);
xlabel('Eb/N0(dB)');
ylabel('BER');
grid on;