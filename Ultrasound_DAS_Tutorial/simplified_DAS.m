%% Delay and Sum beamformer for Focused B-Mode Ultrasound Imaging
clc
clear all
close all

%% Load US-RF data 
load US_DATA_SIMPLE.mat
% US_DATA_SIMPLE = US_DATA_SIMPLE(:,50:end,:);
% size/dimensions of data
N_sl        = size(US_DATA_SIMPLE,1);    % number of scanlines
DepthSample = size(US_DATA_SIMPLE,2);    % samples along axis dimension
N_ele       = size(US_DATA_SIMPLE,3);    % number of elements


%% System parameters
% physical constants
fs = 40e6;               % sampling frequency (Hz)
c = 1540;                % speed of sound in water (m/sec)
pitch = 0.02 * 1e-2;     % cm => m
axial_res = c/fs/2;      % axial resolution c = 2*d/t ==> 2*d*fs (m) (from round trip distance relation)

% physical position information of samples
axial_pos = [0:DepthSample-1]' * axial_res; % index to distance conversion
d_mtx = axial_pos * ones(1,N_ele); % depth information matrix

% probe/elements specs (physical position and size)
scan_view_size = pitch*N_ele; % Lateral scan view size
sc_d = scan_view_size/(N_sl); % Scanline distance/density/resolution/seperation (lateral size of covered by a scanline)
st_sc_x = - scan_view_size*0.5+sc_d*0.5 ; % Start Scanline position (Half of lateral scan view size + half scanline density)

half_ele = -(N_ele-1)*pitch/2; % half of total aperture size
ele_pos = [half_ele:pitch:-half_ele]; % element physical position matrix

% rx data read pointer offset matrix (starting index for each scanline)
rp_mtx = ones(DepthSample,1) * [0:DepthSample:DepthSample*(N_ele-1)];


%% delay-and-sum (DAS) beamformer
beamformer_output  = zeros(DepthSample, N_sl); % focused scanline matrix (uncompressed image)
% Iterating over each scanline   
for sl = 1:N_sl
    disp(['SL : ',num2str(sl)])
    
    % load each scanline at a time
    RF_Ch_data = squeeze(US_DATA_SIMPLE(sl,:,:));

   %%%%%%%%%%%% focusing channel data (apply delays) %%%%%%%%%%%%%%%%%%% 
    sl_pos = st_sc_x + (sl-1)*sc_d; % current scanline position (lateral axis)
    ch_pos = abs(sl_pos - ele_pos); % current active channels position
    ch_pos_mtx = ones(DepthSample,1) * ch_pos;  % channel position matrix

    %% Write your code in this section
    % -------------------------------------%
    tx_d = d_mtx;  % distance from transmitter to reflectivity point (image-pixel/object]
    rx_d = sqrt(ch_pos_mtx.^2 + d_mtx.^2); % distance from reflectivity point (image-pixel/object] to receiver
    
    tx_t = tx_d./c;  % time of flight emitter to point of reflectivity
    rx_t = rx_d./c;  % reflective point to reciever 
    
    tof_rip = tx_t + rx_t; % round trip time of flight
    % -------------------------------------%
    
    % convert delays into index
    read_pointer_rx = round((tof_rip)*fs);  % (tof_rip x fs ) -> index

    % check whether the desired point is in the axial range or not and use edge sample for outliers
    idx = find(read_pointer_rx > DepthSample);
    read_pointer_rx(idx) = DepthSample; 

    idx = find(read_pointer_rx < 1);
    read_pointer_rx(idx) = 1;
    
    % apply delays to RF_Ch_data
    focused_ch_data = RF_Ch_data(read_pointer_rx+rp_mtx)';
    
    %%%%%%%%%%%% beamformer (sum) %%%%%%%%%%%%%%%%%%%
    %% Write your code in this section
    % -------------------------------------%
    beamformer_output(:, sl) = sum(focused_ch_data,1);
%     beamformer_output(:, sl) = mean(focused_ch_data).*((sum(focused_ch_data).^2)./sum(focused_ch_data.^2));
    % -------------------------------------%
    focused_ch_data_mtx(sl,:,:) =focused_ch_data;
    pause(0.05)
    figure(1), subplot(1,2,1), imagesc(log10(abs(focused_ch_data'))), colormap gray, title([int2str(sl),'SCANLINE DATA'])
    figure(1), subplot(1,2,2), imagesc(log10(abs(beamformer_output))), colormap gray, title('B-MODE IMAGE')

end    

% for ch=1:size(focused_ch_data_mtx,2)
%     pause(0.5)
%     figure(1), imagesc(log10(abs(squeeze(focused_ch_data_mtx(:,ch,:))'))), colormap gray, title([int2str(sl),'IMAGE'])
% end
%% Envelope detection
% Envelope detection : in-phase (I) and quadrature (Q) - IQ to image
%% Write your code in this section
% -------------------------------------%
analytic_signal = hilbert(beamformer_output); % calculate analytic function using hilbert transform
RF_env = abs(analytic_signal); % from analytic signal detect envelope signal
% -------------------------------------%

%% Log-compression
dB = 60;                            % dynamic-range in dB
min_dB = 10^(-dB/20);               % zero threshold
norm_data = RF_env./max(RF_env(:)); % scaling to 0 ~ 1 range

%% Write your code in this section
% -------------------------------------%
for i=1:N_sl
    for j=1:DepthSample
        if(norm_data(j,i) < min_dB)
            log_data(j,i) = -dB;  % zero threshold
        else
            log_data(j,i) = ((20)*log10(norm_data(j,i))); % compression
        end
    end
end
% -------------------------------------%

%% visualize image
lateral_tick = 1000*sc_d*[0:N_sl-1];
axial_tick = 1000*axial_res*[0:DepthSample-1];
figure, imagesc(lateral_tick,axial_tick,log_data),colormap gray
c = colorbar; c.Label.String = 'dB'; 
xlabel('Lateral length (mm)')
ylabel('Axial depth (mm)')
% saveas(gcf,'Thyroid_DAS_60dB.bmp','bmp');

%% Calculate Contrast Recovery
% selecting (hypoechoic/anechoic regon) from area in the dimensions ([19,23] to [20.5 23])
cyst_axial_index = round([19 23]./(1000*axial_res));
cyst_lateral_index = round([20.5 23]./(1000*sc_d));
% cyst_lateral_index = round([4 8]./(1000*sc_d));
cyst_mask = zeros(size(log_data));
cyst_mask(cyst_axial_index(1):cyst_axial_index(2),cyst_lateral_index(1):cyst_lateral_index(2))=1;
cyst_region = log_data(find(cyst_mask==1));

% selecting (hyperechoic regon) from area in the dimensions ([19,23] to [14.5 17])
background_axial_index = round([19 23]./(1000*axial_res));
background_lateral_index = round([14.5 17]./(1000*sc_d));
% background_lateral_index = round([9.5 13.5]./(1000*sc_d));
background_mask = zeros(size(log_data));
background_mask(background_axial_index(1):background_axial_index(2),background_lateral_index(1):background_lateral_index(2))=1;
background_region = log_data(find(background_mask==1));

% Calculate CR
meu_cyst = mean(cyst_region(:));
meu_back = mean(background_region(:));

Contrast_recovery = abs(meu_back-meu_cyst)

