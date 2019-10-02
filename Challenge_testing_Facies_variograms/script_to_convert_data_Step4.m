
clear all; close all;
PreName = 'HardTest_';
load([PreName, 'DataSyn1_clip_Ricker_Source_final_30Hz.mat']);
%load('DataSyn1_clip_Ormsby_Source_final_30Hz.mat');
% load('StanfordVI_clip_Source_Ormsby30Hz.mat');
% impedance_time = permute(impedance_time,[2,1,3]);
% wz1_with_multiples = permute(wz1_with_multiples, [2,1,3]);

xlength = size(impedance_time, 1);
ylength = size(impedance_time, 2);
depthlength = size(impedance_time, 3);
nTraces = xlength*ylength;

[r,c,v] = ind2sub(size(impedance_time),find(isnan(impedance_time(:))));
index = min(min(v))-1;

% Sample number where source has maximum amplitude
maxWavelet = 25; % 285 Ormsby 30Hz 90 Ricker 30Hz
dumMult = zeros(index,nTraces); % zeros(index,10000);
% dumMult1 = zeros(index,10000);
dumMult2 = zeros(index,nTraces);
for i = 1:xlength
    clc;disp(num2str(i))
    for j=1:ylength
        % index after which impedance does not exist
        
       dumMult(:,(i-1)*ylength+j) = ...
           squeeze(wz1_with_multiples(i,j,maxWavelet:maxWavelet+index-1));
%        dumMult1(:,(i-1)*10+j) = squeeze(reflect_coeff(i,j,1:index));
%        z0 = RHOBMat(i,j,1)*VpMat(i,j,1);
%        dumMult2(:,(i-1)*100+j) = z0*exp(2*cumsum(reflect_coeff(i,j,:)));
       dumMult2(:,(i-1)*ylength+j) = squeeze(impedance_time(i,j,1:index));
    end
end


wz_with_multiples = dumMult;
% reflect_coeff = dumMult1;
IpTimeVec = dumMult2;
save([PreName, 'DataSyn_Ricker30.mat'],'wz_with_multiples','IpTimeVec')
%save('StanfordVI_Ormsby30.mat','wz_with_multiples','IpTimeVec')


%% Plotting Ip vs waveform

close all;

figure;


sampleNumber = 20;
yyaxis left
plot(IpTimeVec(:,sampleNumber)./1e6);
hold on;
yyaxis right
plot(wz_with_multiples(:,sampleNumber).*3);
title('seismic vs. impedance');
