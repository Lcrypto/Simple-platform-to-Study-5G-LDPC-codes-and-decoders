load('ber_nms_064_054.mat');

figure;

semilogy(list_snr,list_ber(1,:),'bx-');
hold on;
semilogy(list_snr,list_ber(2,:),'r+-');
semilogy(list_snr,list_ber(3,:),'m^-');
semilogy(list_snr,list_ber(4,:),'kv-');
semilogy(list_snr,list_ber(5,:),'gd-');
semilogy(list_snr,list_ber(6,:),'co-');

legend('\alpha=0.64',...
    '\alpha=0.62',...
    '\alpha=0.60',...
    '\alpha=0.58',...
    '\alpha=0.56',...
    '\alpha=0.54'...
    );

axis([-2 -1 10e-5 1]);
xlabel('Eb/N0(dB)');
ylabel('BER');
grid on;