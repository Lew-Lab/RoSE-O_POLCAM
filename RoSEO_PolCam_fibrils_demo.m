clear; close all;
addpath('utils');

%% load raw data

% camera properties
pxSize = 57.5;
offset = 100;
gain = 2;

for i = 1:1e3
    img_n45(:,:,i) = imread('fibrils_NileRed_n45.tif',i);
    img_0(:,:,i) = imread('fibrils_NileRed_0.tif',i);
    img_90(:,:,i) = imread('fibrils_NileRed_90.tif',i);
    img_p45(:,:,i) = imread('fibrils_NileRed_p45.tif',i);
end

img_n45 = flipud(double(img_n45(105:204,170:269,:))-offset)*gain;
img_0 = flipud(double(img_0(105:204,170:269,:))-offset)*gain;
img_90 = flipud(double(img_90(105:204,170:269,:))-offset)*gain;
img_p45 = flipud(double(img_p45(105:204,170:269,:))-offset)*gain;
img = [img_n45,img_0;img_90,img_p45];
disp('raw data loaded');

%% background estimation

b_n45 = GaussWavelet_backgroundFit(mean(img_n45,3));
b_0 = GaussWavelet_backgroundFit(mean(img_0,3));
b_90 = GaussWavelet_backgroundFit(mean(img_90,3));
b_p45 = GaussWavelet_backgroundFit(mean(img_p45,3));

b = [b_n45,b_0;b_90,b_p45];
b = b(1:end-1,1:end-1); % RoSE-O works with images with odd size
b_mean = mean(b(:));
disp('background estimation complete');

%% compute basis images
imgSz = size(img_n45,1);
PSFsz_small = 15;
[FPSF,Bs,B_Ha,B,Bgradx,Bgrady,sumNorm] = computeBasis(imgSz,PSFsz_small,'on');
disp('basis generated');

%% RoSE-O estimation

reg = .25; % regularizer
loc_data = [];
obj.pixelUpsample = 1;
obj.pixelSize = pxSize*2;

% number of workers in parfor
parWorker = feature('numcores');
for groupInd = 1:ceil(size(img_n45,3)/parWorker)
    % RoSE-O works with images with odd size
    if groupInd ~= ceil(size(img_n45,3)/parWorker) || rem(size(img,3),parWorker) == 0
        SMLM_img = img(1:end-1,1:end-1,groupInd*parWorker-parWorker+(1:parWorker));
    else
        SMLM_img = img(1:end-1,1:end-1,groupInd*parWorker-parWorker+(1:min([parWorker,rem(size(img_n45,3),parWorker)])));
    end
    parfor subInd = 1:size(SMLM_img,3)
        [~,~,locDataTmp] = RoSEO_PolCam(obj, SMLM_img(:,:,subInd), b, FPSF, 'regVal', reg);
        if ~isempty(locDataTmp)
            locDataTmp(:,1) = groupInd*parWorker+subInd-parWorker; % record frame number
            for i = 1:size(locDataTmp,1)
                % project to first order moments
                [mux, muy, muz, rotMobil] = secondM2SymmConeWeighted(Bs,B_Ha,sumNorm,locDataTmp(i,5:10),locDataTmp(i,4),b_mean);
                locDataTmp(i,11:14) = [mux, muy, muz, rotMobil]; % -muy to match coordinates
%                 locDataTmp(i,2) = -locDataTmp(i,2);
            end
            loc_data = [loc_data;locDataTmp];
        end
    end
    disp(['frames ',num2str(groupInd*parWorker-parWorker+1),'-',num2str(min([groupInd*parWorker,size(img_n45,3)])),' analyzed']);
    % visualization
    locTmp = loc_data(loc_data(:,1) == loc_data(end,1),[2,3]);
    figure(101);
    subplot(2,2,1);
    imagesc((-imgSz:2:imgSz)*pxSize,(-imgSz:2:imgSz)*pxSize,img_n45(:,:,loc_data(end,1))-b_n45); hold on; 
    plot(locTmp(:,1),locTmp(:,2),'rx'); hold off; axis image;
    subplot(2,2,2);
    imagesc((-imgSz:2:imgSz)*pxSize,(-imgSz:2:imgSz)*pxSize,img_0(:,:,loc_data(end,1))-b_0); hold on;
    plot(locTmp(:,1),locTmp(:,2),'rx'); hold off; axis image;
    subplot(2,2,3);
    imagesc((-imgSz:2:imgSz)*pxSize,(-imgSz:2:imgSz)*pxSize,img_90(:,:,loc_data(end,1))-b_90); hold on;
    plot(locTmp(:,1),locTmp(:,2),'rx'); hold off; axis image;
    subplot(2,2,4);
    imagesc((-imgSz:2:imgSz)*pxSize,(-imgSz:2:imgSz)*pxSize,img_p45(:,:,loc_data(end,1))-b_p45); hold on;
    plot(locTmp(:,1),locTmp(:,2),'rx'); hold off; axis image;
end

save('fibrils_NileRed_flip.mat','loc_data','b_n45','b_0','b_90','b_p45');

%% visualize final results
sideL = pxSize*imgSz;
binSz = 80;
visualizeRoSEest(loc_data,sideL,binSz,'off')

