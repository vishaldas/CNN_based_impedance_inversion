%%
% cd('C:\Users\uriwolln\Desktop\cs229\Exported_properties')
%% run Kennett 2D
close all;
cd ('/Volumes/GoogleDrive/My Drive/Research/CNN_seismic_inversion/Code_Stuff/Facies_variograms/ABC_dataset');
variogramRanges = linspace(10, 60, 21);
% variogramRanges = [3 5 8 70 80 90 100];
% variogramRanges = [40];
PreName = 'Train_';

for i=1:length(variogramRanges)
    i
    clearvars -except i zClayMat zPhiMat zFaciesMat variogramRanges PreName
    variogramRange = variogramRanges(i);
    
    xlength = 10;
    ylength = 10;
    zlength = 200;
    
    zFaciesMat_temp = zeros(xlength,ylength,zlength);
    zClayMat_temp = zeros(xlength,ylength,zlength);
    zPhiMat_temp = zeros(xlength,ylength,zlength);

    slice_n_cells = xlength*ylength;
    fid = fopen(['facies_clay_porosity_', num2str(variogramRange),'.txt'],'r');
    tline = fgetl(fid);
    for whichSlice = 1:zlength %+800
        %whichSlice
        zFacies = zeros(slice_n_cells,1);
        zClay = zeros(slice_n_cells,1);
        zPhi = zeros(slice_n_cells,1);
        for j=1:slice_n_cells
            tline = fgetl(fid);
            dum = str2num(tline);
            zFacies(j) = dum(1);
            zClay(j) = dum(2);
            zPhi(j) = dum(3);
        end
        I = (reshape(zFacies, [xlength  ylength]));
        zFaciesMat_temp(:,:,whichSlice) = I;
        
        I = (reshape(zClay, [xlength  ylength]));
        zClayMat_temp(:,:,whichSlice) = I;
        
        I = (reshape(zPhi, [xlength  ylength]));
        zPhiMat_temp(:,:,whichSlice) = I;
        
        %     I = (reshape(zFacies, [xlength  ylength]));
        %     I = I';
        %     figure;imagesc(flipud(I))
    end
      
%     figure
%     subplot(1,3,1)
%     plot(squeeze(zFaciesMat_temp(1,1,:)), 1:200)
%     xlim([-2, 2])
%     title(['facies with variogram ', num2str(variogramRange)]);
%     set(gca,'Ydir','reverse')
%     
%     subplot(1,3,2)
%     plot(squeeze(zPhiMat_temp(1,1,:)), 1:200)
%     xlim([0, 1])
%     title('phi');
%     set(gca,'Ydir','reverse')
%     
%     subplot(1,3,3)
%     plot(squeeze(zClayMat_temp(1,1,:)), 1:200)
%     xlim([0, 1])
%     title('clay');
%     set(gca,'Ydir','reverse')
    
    if(i==1)
        zClayMat = zClayMat_temp;
        zPhiMat = zPhiMat_temp;
        zFaciesMat = zFaciesMat_temp;
    else
        zClayMat = [zClayMat zClayMat_temp];
        zPhiMat = [zPhiMat zPhiMat_temp];
        zFaciesMat = [zFaciesMat zFaciesMat_temp];
    end
end

    save([PreName 'ClayMat.mat'],'zClayMat')
    save([PreName 'PhiMat.mat'],'zPhiMat')
    save([PreName 'FaciesMat.mat'],'zFaciesMat')
