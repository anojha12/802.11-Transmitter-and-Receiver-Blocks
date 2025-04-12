
function [decoded_data]= MyOfdmReceiver(data);

OFDM_TX; 
 
%% Rx processing params
rx_data = raw_rx_data;         
LTS_CORR_THRESH = 0.5;  
STS_CORR_THRESH = 0.5;

%% Packet Detection using STS: 
%Auto correlation of received signal with STS to detect the packet:


for i =  1:length(rx_data)-32;
    window1= rx_data(i:i+15);
    window2= rx_data(i+16:i+31);
    output_auto_sts(i)= (dot(window1,(window2)'))/(norm(window1)*norm(window2) + 1e-6);
end

fprintf("First peak detected using sts_auto_correlation: %d\n", find(output_auto_sts > STS_CORR_THRESH, 1, 'first'));
figure;
plot(abs(output_auto_sts));  
hold on;
[peaks, locs] = findpeaks(abs(output_auto_sts), 'MinPeakHeight', STS_CORR_THRESH);
plot(locs, peaks, 'go', 'MarkerFaceColor', 'r');
title('Auto-Correlation with STS');
grid on;
xlim([0, 1000]); 

%Cross correlation of received signal with STS to detect the packet:

for i =  1:length(rx_data)-16;
    window1= rx_data(i:i+15);
    output_cross_sts(i)= (dot(window1,sts_t))/(norm(window1)*norm(sts_t) + 1e-6);
end


fprintf("First peak detected using sts_cross_correlation: %d\n", find(output_cross_sts > STS_CORR_THRESH, 1, 'first'));
figure;
plot(abs(output_cross_sts));  
hold on;
[peaks, locs] = findpeaks(abs(output_cross_sts), 'MinPeakHeight', STS_CORR_THRESH);
plot(locs, peaks, 'go', 'MarkerFaceColor', 'r');
title('Cross-Correlation with STS');
grid on;
xlim([0, 1000]);  

%Cross correlation of received signal with LTS to detect the packet:

for i =  1:length(rx_data)-64;
    window1= rx_data(i:i+63);
    output_cross_lts(i)= (dot(window1,lts_t))/(norm(window1)*norm(lts_t) + 1e-6);
end


fprintf("First peak detected using lts_cross_correlation: %d\n", find(output_cross_lts > LTS_CORR_THRESH, 1, 'first'));
figure;
plot(abs(output_cross_lts));  
hold on;
[peaks, locs] = findpeaks(abs(output_cross_lts), 'MinPeakHeight', LTS_CORR_THRESH);
plot(locs, peaks, 'go', 'MarkerFaceColor', 'r');
title('Cross-Correlation with LTS');
grid on;
xlim([0, 1000]);

%% CFO estimation and correction
% Use two copies of LTS for cross-correlation (Reference: Thesis)

sts_pkt_start= find(output_cross_sts > STS_CORR_THRESH, 1, 'first');
sts_pkt_end= sts_pkt_start + 30*16;
lts_pkt_start= sts_pkt_end +32;
data_pkt_start= lts_pkt_start+128;

fprintf("sts_pkt_start %d\n",sts_pkt_start);
fprintf("sts_pkt_end %d\n",sts_pkt_end);
fprintf("lts_pkt_start %d\n",lts_pkt_start);
fprintf("lts_pkt_end %d\n",data_pkt_start);

% Output: Packet with each value multiplied by CFO correction factor
%Coarse CFO
for i = sts_pkt_start:sts_pkt_end-32
    window1= rx_data(i:i+15);
    window2= rx_data(i+16:i+31);
    coarse_cfo = sum(window2.*conj(window1));
end

coarse_cfo= abs(angle(coarse_cfo)/(-2*pi*16));
fprintf("Coarse CFO is: %d\n",coarse_cfo);

%Residual CFO
window1= rx_data(lts_pkt_start:lts_pkt_start+63);
window2= rx_data(lts_pkt_start+64:data_pkt_start-1);
residual_cfo= sum(window2 .* window1);
residual_cfo= abs(angle(residual_cfo)/(-2*pi*64));

fprintf("Residual CFO is: %d\n",residual_cfo);

total_cfo= coarse_cfo+residual_cfo;
fprintf("Total CFO is: %d\n",total_cfo);

%CFO Correction:
rx_data= rx_data.*exp(-j*2*pi*total_cfo*length(rx_data));

%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help
% Output: CP free payload matrix of size (N_SC * N_OFDM_SYMS)
rx_syms_data= rx_data(data_pkt_start:data_pkt_start+(N_SC+CP_LEN)*N_OFDM_SYMS-1);
rx_mat= reshape(rx_syms_data,(N_SC+CP_LEN),[]);
nocp_rx= rx_mat(CP_LEN+1:end,:);


%% FFT
% Refer to IFFT perfomed at TX
% Output: Symbol matrix in frequency domain of same size
rx = fft(nocp_rx,N_SC,1);

%% Channel estimation and correction
% Use the two copies of LTS and find channel estimate (Reference: Thesis)
% Convert channel estimate to matrix form and equlaize the above matrix
% Output : Symbol equalized matrix in frequency domain of same size
lts_f1= fft(rx_data(lts_pkt_start:lts_pkt_start+63),64);
lts_f2= fft(rx_data(lts_pkt_start+64:data_pkt_start-1),64);
lts_f1(lts_f1==0)= 1;
lts_f2(lts_f2==0)= 1;
lts_f(lts_f==0)= 1;

hest1= lts_f1 ./ (lts_f);
hest2= lts_f2 ./ (lts_f);
hest= (hest1+hest2)/2;
hest(hest==0)= 1;

figure;
plot(abs(hest))
xlabel("Frequency")
ylabel("Estimated Channel")
title("Estimated Channel")


rx= rx./conj(hest');

%% Advanced topics: 
%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.

% Output: Symbol equalized matrix with pilot phase correction applied


for i= 1:N_OFDM_SYMS
    rx_tray= rx(:, i);
    pilots_rx= rx_tray(SC_IND_PILOTS);
    pilots_rx_reorder = [pilots_rx(1) pilots_rx(2) pilots_rx(3) pilots_rx(4)]';
    pilot_phase_shift = angle(pilots_rx_reorder ./ pilots);
    residual_cfo = mean(pilot_phase_shift);
    rx(:, i) = rx(:, i) * exp(1j * residual_cfo);
end

%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error

% Output: Symbol equalized matrix with pilot phase correction applied
% Remove pilots and flatten the matrix to a vector rx_syms
rxdata= rx(SC_IND_DATA,:); 
rx_syms= reshape(rxdata, 1, []);

%% Demodulation

figure(4);
scatter(real(rx_syms), imag(rx_syms),'filled');
title(' Signal Space of received bits');
xlabel('I'); ylabel('Q');

decoded_data = [];
% FEC decoder
Demap_out = demapper(rx_syms,MOD_ORDER,1);

% viterbi decoder
decoded_data_func = vitdec(Demap_out,trel,7,'trunc','hard');
decoded_data = [decoded_data;decoded_data_func];
% decoded_data is the final output corresponding to tx_data, which can be used
% to calculate BER

%BER Calculation:
BER=  sum(decoded_data ~= tx_data);
fprintf("BER is: %d\n",(BER/length(tx_data)*100));
