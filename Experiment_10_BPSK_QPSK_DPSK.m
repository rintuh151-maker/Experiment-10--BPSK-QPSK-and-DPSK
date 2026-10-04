clc;
clear;
close all;

%% EXPERIMENT 10 - BPSK, QPSK AND DPSK
% Normalized BPSK, Gray-coded QPSK and DBPSK
% Phase offsets: 10, 30, 60 and 90 degrees

rng(10);

N = 100000;
EbN0_dB = 0:2:12;
phase_offsets = [10 30 60 90];

%% ---------------------------------------------------------
% 1. BPSK
%% ---------------------------------------------------------

bits_bpsk = randi([0 1],1,N);

% BPSK mapping: 0 -> +1, 1 -> -1
tx_bpsk = 1 - 2*bits_bpsk;

%% ---------------------------------------------------------
% 2. Gray-coded QPSK
%% ---------------------------------------------------------

Nq = 100000;
bits_qpsk = randi([0 1],1,2*Nq);

% Gray mapping:
% 00 -> ( +1 + j )/sqrt(2)
% 01 -> ( -1 + j )/sqrt(2)
% 11 -> ( -1 - j )/sqrt(2)
% 10 -> ( +1 - j )/sqrt(2)

b1 = bits_qpsk(1:2:end);
b2 = bits_qpsk(2:2:end);

I = 1 - 2*b2;
Q = 1 - 2*b1;

tx_qpsk = (I + 1j*Q)/sqrt(2);

%% ---------------------------------------------------------
% 3. Differential encoding and DBPSK
%% ---------------------------------------------------------

bits_dpsk = randi([0 1],1,N);

% Differential encoding:
% bit 0 -> no phase change
% bit 1 -> 180 degree phase change

differential_bits = zeros(1,N);
differential_bits(1) = bits_dpsk(1);

for k = 2:N
    differential_bits(k) = xor(differential_bits(k-1),bits_dpsk(k));
end

tx_dpsk = 1 - 2*differential_bits;

%% ---------------------------------------------------------
% 4. Mandatory manual validation
%% ---------------------------------------------------------

disp('---------------------------------------------');
disp('MANDATORY QPSK MAPPING VALIDATION');
disp('---------------------------------------------');

manual_bits = [0 0; 0 1; 1 1; 1 0];

for k = 1:4
    b1m = manual_bits(k,1);
    b2m = manual_bits(k,2);

    Im = 1 - 2*b2m;
    Qm = 1 - 2*b1m;

    symbol = (Im + 1j*Qm)/sqrt(2);

    fprintf('%d%d -> I = %+d, Q = %+d -> Symbol = %.3f %+.3fj\n', ...
        b1m,b2m,Im,Qm,real(symbol),imag(symbol));
end

%% ---------------------------------------------------------
% 5. Ideal / noisy constellations
%% ---------------------------------------------------------

SNR_const = 8;

noise_qpsk = sqrt(1/(2*10^(SNR_const/10))) * ...
    (randn(1,Nq) + 1j*randn(1,Nq));

rx_qpsk_noisy = tx_qpsk + noise_qpsk;

figure;
plot(real(tx_qpsk(1:3000)),imag(tx_qpsk(1:3000)),'.');
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('QPSK Ideal Constellation');

figure;
plot(real(rx_qpsk_noisy(1:3000)),imag(rx_qpsk_noisy(1:3000)),'.');
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('QPSK Noisy Constellation');

%% ---------------------------------------------------------
% 6. QPSK mapping visualization
%% ---------------------------------------------------------

figure;
hold on;

map_symbols = [
    (1+1j)/sqrt(2)
    (-1+1j)/sqrt(2)
    (-1-1j)/sqrt(2)
    (1-1j)/sqrt(2)
];

plot(real(map_symbols),imag(map_symbols),'o','MarkerSize',8);

text(real(map_symbols(1))+0.05,imag(map_symbols(1)),'00');
text(real(map_symbols(2))+0.05,imag(map_symbols(2)),'01');
text(real(map_symbols(3))+0.05,imag(map_symbols(3)),'11');
text(real(map_symbols(4))+0.05,imag(map_symbols(4)),'10');

grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Gray-Coded QPSK Mapping');
hold off;

%% ---------------------------------------------------------
% 7. Phase trajectory
%% ---------------------------------------------------------

phase_data = unwrap(angle(tx_dpsk(1:200)));

figure;
plot(1:200,phase_data,'-');
grid on;
xlabel('Symbol Number');
ylabel('Phase (radians)');
title('DBPSK Phase Trajectory');

%% ---------------------------------------------------------
% 8. BER - coherent BPSK with phase offsets
%% ---------------------------------------------------------

BER_BPSK = zeros(length(phase_offsets),length(EbN0_dB));

for p = 1:length(phase_offsets)

    theta = phase_offsets(p)*pi/180;

    for s = 1:length(EbN0_dB)

        noise = sqrt(1/(2*10^(EbN0_dB(s)/10))) * ...
            randn(1,N);

        rx = tx_bpsk*exp(1j*theta) + noise;

        detected = real(rx) < 0;

        BER_BPSK(p,s) = mean(detected ~= bits_bpsk);

    end
end

%% ---------------------------------------------------------
% 9. BER - coherent Gray-coded QPSK
%% ---------------------------------------------------------

BER_QPSK = zeros(1,length(EbN0_dB));

for s = 1:length(EbN0_dB)

    noise = sqrt(1/(2*10^(EbN0_dB(s)/10))) * ...
        (randn(1,Nq) + 1j*randn(1,Nq));

    rx = tx_qpsk + noise;

    detected_b1 = imag(rx) < 0;
    detected_b2 = real(rx) < 0;

    detected_bits = zeros(1,2*Nq);

    detected_bits(1:2:end) = detected_b1;
    detected_bits(2:2:end) = detected_b2;

    BER_QPSK(s) = mean(detected_bits ~= bits_qpsk);

end

%% ---------------------------------------------------------
% 10. BER - DBPSK differential detection
%% ---------------------------------------------------------

BER_DPSK = zeros(1,length(EbN0_dB));

for s = 1:length(EbN0_dB)

    noise = sqrt(1/(2*10^(EbN0_dB(s)/10))) * ...
        (randn(1,N) + 1j*randn(1,N));

    rx = tx_dpsk + noise;

    differential_product = rx(2:end).*conj(rx(1:end-1));

    detected = real(differential_product) < 0;

    BER_DPSK(s) = mean(detected ~= bits_dpsk(2:end));

end

%% ---------------------------------------------------------
% 11. Phase offset comparison
%% ---------------------------------------------------------

BER_offset = zeros(length(phase_offsets),length(EbN0_dB));

for p = 1:length(phase_offsets)

    theta = phase_offsets(p)*pi/180;

    for s = 1:length(EbN0_dB)

        noise = sqrt(1/(2*10^(EbN0_dB(s)/10))) * ...
            (randn(1,N) + 1j*randn(1,N));

        rx = tx_dpsk*exp(1j*theta) + noise;

        differential_product = rx(2:end).*conj(rx(1:end-1));

        detected = real(differential_product) < 0;

        BER_offset(p,s) = mean(detected ~= bits_dpsk(2:end));

    end
end

%% ---------------------------------------------------------
% 12. BER comparison
%% ---------------------------------------------------------

figure;
semilogy(EbN0_dB,BER_BPSK(1,:),'-o');
hold on;
semilogy(EbN0_dB,BER_QPSK,'-s');
semilogy(EbN0_dB,BER_DPSK,'-^');

grid on;
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Comparison');
legend('BPSK','Gray-coded QPSK','DBPSK','Location','southwest');
hold off;

%% ---------------------------------------------------------
% 13. Expected physical effects
%% ---------------------------------------------------------

disp(' ');
disp('---------------------------------------------');
disp('OBSERVATION AND INTERPRETATION');
disp('---------------------------------------------');

disp('Expected effect: Increasing phase offset rotates the received constellation.');
disp('Expected effect: Coherent detection becomes more sensitive to phase error.');
disp('Expected effect: Differential detection removes the need for absolute carrier phase.');
disp('Expected effect: Increasing Eb/N0 should reduce BER.');

disp(' ');
disp('Phase offsets tested:');
disp('10, 30, 60 and 90 degrees');

disp(' ');
disp('Experiment 10 completed.');