load('ber_oms_050_060.mat');

figure;

semilogy(list_snr,list_ber(1,:),'bx-');
hold on;
semilogy(list_snr,list_ber(2,:),'r+-');
semilogy(list_snr,list_ber(3,:),'m^-');
semilogy(list_snr,list_ber(4,:),'kv-');
semilogy(list_snr,list_ber(5,:),'gd-');
semilogy(list_snr,list_ber(6,:),'co-');

legend('\beta=0.50',...
    '\beta=0.52',...
    '\beta=0.54',...
    '\beta=0.56',...
    '\beta=0.58',...
    '\beta=0.60'...
    );

axis([-2 -1 10e-5 1]);
xlabel('Eb/N0(dB)');
ylabel('BER');
grid on;