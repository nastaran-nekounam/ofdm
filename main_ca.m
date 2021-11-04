% clear all
rng(1);
% set parameters
bit_num=1e7;
ifft_size = 1024;
nc = 400;
symb_size = 2;
clipping = 3;
%SNR_dB = 40;
guard_time = ifft_size/4;
spf = ceil(2^13/nc);
symb_period = ifft_size + guard_time;
head_len = symb_period*8;
envelope = ceil(symb_period/256)+1;
spacing = 0;
while (nc*spacing) <= (ifft_size/2 - 2)
    spacing = spacing + 1;
end
spacing = spacing - 1;
midFreq = ifft_size/4;
FC = midFreq - round((nc-1)*spacing/2);
LC = midFreq + floor((nc-1)*spacing/2);
carriers = [FC:spacing:LC] + 1;
conj_carriers = ifft_size - carriers + 2;
save('ofdm_param');


%*************************sender***************************%
%generate random data
rnd_data_ = randi([0 1],1,bit_num);
save('err_calc.mat', 'rnd_data_');
%qpsk modulation
rnd_data_ = bi2de(reshape(rnd_data_,2,[])');
f = 0.25;
header = sin(0:f*2*pi:f*2*pi*(head_len-1));
f=f/(pi*2/3);
header = header+sin(0:f*2*pi:f*2*pi*(head_len-1));
frame_guard = zeros(1, symb_period);
transmit_data = [];
symb_per_carrier = ceil(length(rnd_data_)/nc);

% frame divider
while ~isempty(rnd_data_)
    
    frame_len = min(spf*nc,length(rnd_data_));
    frame_data = rnd_data_(1:frame_len);
    rnd_data_ = rnd_data_((frame_len+1):(length(rnd_data_)));
    
    carrier_symb_count = ceil(length(frame_data)/nc);
    if length(frame_data)/nc ~= carrier_symb_count,
        padding = zeros(1, carrier_symb_count*nc);
        padding(1:length(frame_data)) = frame_data;
        frame_data = padding;
    end
    % serial to parellel
    data_tx_matrix = reshape(frame_data, nc, carrier_symb_count)';
    %DPSK modulation
    carrier_symb_count = size(data_tx_matrix,1) + 1;
    diff_ref = round(rand(1, nc)*(2^symb_size)+0.5);
    data_tx_matrix = [diff_ref; data_tx_matrix];
    for k=2:size(data_tx_matrix,1)
        data_tx_matrix(k,:) = ...
            rem(data_tx_matrix(k,:)+data_tx_matrix(k-1,:), 2^symb_size);
    end
    [X,Y] = pol2cart(data_tx_matrix*(2*pi/(2^symb_size)),ones(size(data_tx_matrix)));
    complex_matrix = X + 1i*Y;
    %IFFT bins allocatation
    spectrum_tx = zeros(carrier_symb_count, ifft_size);
    spectrum_tx(:,carriers) = complex_matrix;
    spectrum_tx(:,conj_carriers) = conj(complex_matrix);
    %IFFT
    modulated_signal = real(ifft(spectrum_tx'))';
    %cp addition
    end_symb = size(modulated_signal, 2); % end of a symbol period without guard
    modulated_signal = [modulated_signal(:,(end_symb-guard_time+1):end_symb) modulated_signal];
    % parellel to serial
    modulated_signal = modulated_signal';
    modulated_signal = reshape(modulated_signal, 1, size(modulated_signal,1)*size(modulated_signal,2));
    %Cascade frames
    transmit_data = [transmit_data frame_guard modulated_signal];
    frame_power = var(modulated_signal);
end

power = frame_power;
transmit_data = [power*header transmit_data frame_guard power*header];


%*************************channel***************************%
% signal clipping
clipped_peak = (10^(0-(clipping/20)))*max(abs(transmit_data));
transmit_data(find(abs(transmit_data)>=clipped_peak))...
    = clipped_peak.*transmit_data(find(abs(transmit_data)>=clipped_peak))...
    ./abs(transmit_data(find(abs(transmit_data)>=clipped_peak)));

%Rayleigh channel

h = raylrnd(1,1,length(transmit_data));
% transmit_data = ifft(fft(h).*fft(transmit_data));

% channel noise

SNR_linear = 10^(SNR_dB/10);
noise_factor = sqrt(1/SNR_linear);
noise = randn(1,length(transmit_data)) * noise_factor;
recieved_data = transmit_data + noise;

save('received.mat', 'recieved_data','h','noise_factor');

%*************************reciever***************************%

clear all;

load('ofdm_param');
load('received.mat');
recieved_data = recieved_data.';
s_idx = 1;
e_idx = length(recieved_data);
data = [];
phase = [];
is_last = 0;
H = fft(h);

    
unpad = 0;
if rem(bit_num, nc)~=0
    unpad = nc - rem(bit_num, nc);
end
num_frame=ceil((bit_num)*(1/symb_size)/(spf*nc));
for k = 1:num_frame
    % frames detection
    if k==1
        frame_start = s_idx + symb_period + head_len;
        frame_end = min(e_idx,frame_start - 1 + symb_period*(spf+1));
    else
        frame_start = frame_end + 1 + symb_period;
        frame_end = min(e_idx,frame_start - 1 + symb_period*(spf+1));
    end
    
    if k==num_frame
        is_last = 1;
    end
    
    time_wave = recieved_data(frame_start:frame_end);
    
    % serial to parallel
    symb_rx_matrix = reshape(time_wave(1:(symb_period*floor(length(time_wave)/symb_period))),symb_period, floor(length(time_wave)/symb_period));
    %cp removal
    symb_rx_matrix = symb_rx_matrix(guard_time+1:symb_period,:);
    %fft
    rx_spectrum_matrix = fft(symb_rx_matrix)';
    %extract carriers from FFT bins
    rx_spectrum_matrix = rx_spectrum_matrix(:,carriers);
    
    %mmse
%     
%     H_F = H(frame_start:frame_end);
%     H_S = reshape(H_F(1:(symb_period*floor(length(time_wave)/symb_period))),symb_period, floor(length(time_wave)/symb_period));
%     H_S = H_S(guard_time+1:symb_period,:)';
%     H_S = H_S(:,carriers);
%     rx_spectrum_matrix = rx_spectrum_matrix.*conj(H_S)./(H_S.*conj(H_S)+ noise_factor);
%     
%    

    %DPSK demode
    rx_phase = angle(rx_spectrum_matrix)*(180/pi);
    rx_phase = rem((rx_phase+360), 360);
    decoded_phase = diff(rx_phase);
    decoded_phase = rem((decoded_phase+360), 360);
    % parellel to serial
   
    
    decoded_phase = reshape(decoded_phase',1, size(decoded_phase,1)*size(decoded_phase,2));
    decoded_phase = rem((decoded_phase+360), 360);
    base_phase = 360/(2^symb_size);
    data_rx = floor(rem((decoded_phase/base_phase+0.5),(2^symb_size))); 
    if is_last==1
        data_rx = data_rx(1:(length(data_rx)-unpad));
    end
    %cascade frames
    data = [data data_rx];
end
data_rx = data;
% qpsk demodulation
data_out = reshape(de2bi(data_rx)',1,[]);
missed_num =0;
if length(data_out)>(bit_num)
    data_out = data_out(1:bit_num);
elseif length(data_out)<(bit_num)
    missed_num =bit_num-length(data_out);
    fprintf('number of missed = %d\n',missed_num)
end
%*************************err calc***************************%
load('err_calc.mat');
errors = find(rnd_data_(1:end-missed_num)~=data_out);
fprintf('number of errors = %d\n',length(errors))
BER = length(errors)/length(data_out);
fprintf('BER = %f%%\n',BER)
