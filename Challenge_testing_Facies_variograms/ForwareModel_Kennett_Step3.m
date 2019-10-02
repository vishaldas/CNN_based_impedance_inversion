%% Kennett 1D
clear all;
addpath(genpath('/Volumes/GoogleDrive/My Drive/Research/CNN_seismic_inversion/Codes_ML/SeisLab_10.0301'));

% From RP model
pre = 'HardTest_';
load([pre 'DataSyn1_time_clip.mat']);

diffDepth = 1;

xlength = size(VpMat,1);
ylength = size(VpMat,2);
depthlength = size(VpMat,3);
% 
% 
% fdom = [20 30 40];
fdom = 30; % Central frequnecy

for k=1:length(fdom)
%     Fix dt based on what you think
%     dt_nyquist = 1/(2*fdom(k));
%     dt = 0.0351*dt_nyquist;
    dt = 4.3875e-4;
%     N_sample = round(dt_nyquist./dt);
    time = 0:dt:4; % Simulating till 4 secs to prevent wrap-around
    
    if (rem(length(time),2) ~= 0)
        time = time(1:end-1);
    end
    
    % wvlt = sourcewvlt;
    % Ormsby wavelet with freq parameters 4,10,50,70
    wav=s_create_wavelet({'type','min-phase'},{'frequencies',4,10,50,70},{'step',dt*1e3},{'wlength',248});
    wvlt = wav.traces;
    % Padding the end of source with zeros
    wvlt(end:length(time)) = 0;
    wvlt = wvlt.';
    
    n=length(wvlt);
%     om=(2*pi/(n*dt))*[0:1:(n/2 -1)];
    
%     freq=om/(2*pi);
    
    wz_with_multiples = zeros(xlength, ylength, length(time));
%     tz_with_multiples = zeros(100, 100, length(freq));
    wz_without_multiples = zeros(xlength, ylength, length(time));
%     tz_without_multiples = zeros(100, 100, length(freq));
%     refl_down = zeros(100,100,length(time));
    impedance_time = zeros(xlength,ylength,length(time));
    
    % 100 = number of traces in i
    for i=1:xlength
        disp(['Iteration num = ', num2str(i)]);
        for j=1:ylength
            % LYR = [Velocity (m/sec) Density(kg/m^3) Thickness(m)]
            LYR = [squeeze(VpMat(i,j,:))*1000, squeeze(RHOBMat(i,j,:))*1000 , ones(length(VpMat(i,j,:)),1)*diffDepth];
            [wz_with_multiples(i,j,:),~,~] = kennet(LYR,wvlt,dt,2,0,0);
%             [wz_without_multiples(i,j,:),~,tz_without_multiples(i,j,:)] = kennet(LYR,wvlt,dt,0,0,0);
            impedance = LYR(:,1) .* LYR(:,2);
            tt = 2*cumsum(diffDepth./LYR(1:end-1,1));
            tt = [0; tt];
            impedance_interp = interp1(tt, impedance, time, 'next');
            impedance_time(i,j,:) = impedance_interp;
            % Reflection coefficients calculations
%             refl = 0.5*diff(log(impedance_interp));
%             refl = [0 refl];
%             refl(isnan(refl)) = 0;
%             refl_down(i,j,:) = refl;
        end
    end
    
    
    
    % % Saving the full data
    % timeVec = time;
    % depthVec =  [1:length(squeeze(VpMat(1,1,:)))]'*diffDepth;
    % freqVec = freq;
    %
    % save('DataSyn1_full_without_multiples.mat','VpMat','RHOBMat','wz', 'timeVec', 'tz', 'freq', 'depthVec');
    
    
    
    
    % cutting time to 0.4 sec and frequency to 250 HZ
    
    ind1 = find(time >= 0.4, 1);
%     ind2 = find(freq >= 250, 1);
    timeVec1 = time(1:ind1);
    depthVec =  [1:length(squeeze(VpMat(1,1,:)))]'*diffDepth;
%     freq1 = freq(1:ind2);
    wz1 = wz_with_multiples(:,:, 1:ind1);
%     tz1 = tz_with_multiples(:,:, 1:ind2);
    
    wz11 = wz_without_multiples(:,:, 1:ind1);
%     tz11 = tz_without_multiples(:,:, 1:ind2);
%     refl_down1 = refl_down(:,:,1:ind1);
    impedance_time1 = impedance_time(:,:,1:ind1);
    
    wz1_with_multiples = wz1;
%     tz1_with_multiples = tz1;
    
    wz1_without_multiples = wz11;
%     tz1_without_multiples = tz11;
    
%     reflect_coeff = refl_down1;
    impedance_time = impedance_time1;
    
%     load('DataSyn1_clip_Ricker_Source_30Hz.mat', 'wz1_with_multiples', 'wz1_without_multiples');
    
    filename = [pre 'DataSyn1_clip_Ricker_Source_final_' num2str(fdom(k)) 'Hz.mat'];
%     
%     save(filename,'VpMat','RHOBMat','wz1_with_multiples', 'wz1_without_multiples','reflect_coeff','timeVec1', ...
%         'tz1_with_multiples', 'tz1_without_multiples', 'freq1', 'depthVec');
    
    save(filename,'impedance_time','wz1_with_multiples', 'wz1_without_multiples','timeVec1', 'depthVec');
    
       
end
