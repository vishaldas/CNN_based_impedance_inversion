% cd('C:\Users\uriwolln\Desktop\cs229\Exported_properties')
%
clear all;
SwThinkInLayer = 0.23;
Kg = 0.06; 
Kw = 2.61;
RHOw = 1.01;
RHOg = 0.18; 
diffP = 16.5; 
RHOs = 2.65;
PhiC = 0.4; 
Coordination = 14; 
Fudge = 1;

PreName = 'Train_';

load([PreName 'FaciesMat.mat'])
load([PreName 'ClayMat.mat'])
load([PreName 'PhiMat.mat'])

xlength = size(zClayMat,1);
ylength = size(zClayMat,2);
depthlength = size(zClayMat,3);

% fdom = 30;
% dt_nyquist = 1/(2*fdom);
% dt = 0.0351*dt_nyquist;
% N_sample = round(dt_nyquist./dt);
% time = 0:dt:8;
% 
% wvlt = sourcewvlt;
% wvlt(end:length(time)) = 0;
% 
% n=length(wvlt);
% om=(2*pi/(n*dt))*[0:1:(n/2 -1)];
% 
% freq=om/(2*pi);

% wz = zeros(xlength, ylength, length(time));
VpMat = zeros(xlength,ylength,depthlength);
RHOBMat = zeros(xlength,ylength,depthlength);

% tz = zeros(xlength, ylength, length(freq));


for i =1:xlength
    disp([num2str(i) '/' num2str(xlength)])
    for j=1:ylength
        
%         load('FaciesMat.mat')
        Facies = squeeze(zFaciesMat(i,j,:));
%         clear zFaciesMat
%         load('PhiMat.mat')
        Phi = squeeze(zPhiMat(i,j,:));
%         clear zPhiMat
%         load('ClayMat.mat')
        Clay = squeeze(zClayMat(i,j,:));
%         clear zClayMat
        
        
        diffDepth = 0.1524;
        depth = [1:length(Clay)]'*diffDepth;
       
        Facies(1) = 0;
        Facies(2) = 1;
        Facies(3) = 0;
        xx = diff(Facies);
        [indx1] = find(xx == 1);
        [indx_min1] = find(xx == -1);
        m = length(indx1);
        n = length(indx_min1);
        l = min([m n]);
        
        if indx1(1) < indx_min1(1)  %%% first layer is shale
            x = indx_min1(1:l) - indx1(1:l) ;
            indx = find(x == max(x));
            try
                indxReservoir = (indx1(indx)+1):(indx1(indx+1)+1);
            catch
                indxReservoir = (indx1(indx)+1):length(Sw);
            end
        else
            x = indx1(2:l) -indx_min1(1:(l-1)) ;
            indx = find(x == max(x));
            indxReservoir = (indx1(indx)+1):(indx1(indx+1)+1);
        end
        
        try
                endGas = round(rand(1)*length((indx1(indx)+1):(indx1(indx+1)+1)) + (indx1(indx)+1));
        
            catch
                endGas = round(rand(1)*length((indx1(indx)+1):(depthlength-1)))+ (indx1(indx)+1);
        
        end
        
        Sw = ones(size(Clay));
        Sw((indx1(indx)+1):endGas) = normrnd(SwThinkInLayer,0.03,length((indx1(indx)+1):endGas),1);
        
        RHOf = Sw*RHOw + (1 - Sw) * RHOg;
        Kf = 1./(Sw./Kw + (1 - Sw)./Kg);
        
        %     Kf = (SwDum/Kw + (1 - SwDum)/Kg)^(-1);
        
        [Vp,Vs,~,~,RHOB]=SoftSandNew(Kf,RHOf,Phi,(1 - Clay)*0.8,Clay, (1 - Clay)*0.2,0,0,diffP, PhiC, Coordination, Fudge,{'c',{10.5 3.5 2.65}});
        if(sum(isnan(Vp))>0)
            Vp = repnan(Vp,'next'); 
            Vp = repnan(Vp,'previous');
        end
        
%         
%         figure;
%         subplot 151
%         plot(Facies,depth); axis ij
%         subplot 152
%         plot(Phi,depth); axis ij
%         subplot 153
%         plot(Clay,depth); axis ij
%         subplot 154
%         plot(Sw,depth); axis ij
%         subplot 155
%         plot(Vp,depth); axis ij ; hold on;  plot(Vs,depth);
        
        LYR = [Vp.*1000, RHOB.*1000, ones(size(Vp))*diffDepth];
       
        
        VpMat(i,j,:) = Vp;
        RHOBMat(i,j,:) = RHOB;
        
    end
end




parfor i=1:1
    for j=1:1
        LYR = [squeeze(VpMat(i,j,:))*1000, squeeze(RHOBMat(i,j,:))*1000 , ones(length(VpMat(1,1,:)),1)*diffDepth];
%         [wz(i,j,:),~,tz(i,j,:)] = kennet(LYR,wvlt,dt,2,0,0);
    end
end

% timeVec = time;
depthVec =  [1:length(squeeze(VpMat(1,1,:)))]'*diffDepth;
% freqVec = freq;

save([PreName 'DataSyn1_time_clip.mat'],'VpMat','RHOBMat', 'depthVec');


