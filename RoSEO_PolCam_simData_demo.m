clear; close all;
addpath('utils');

%% compute basis images
imgSz = 25;
PSFsz_small = 15; % for projection to 1st orientational moments
[FPSF,Bs,B_Ha,B,Bgradx,Bgrady,sumNorm] = computeBasis(imgSz,PSFsz_small,'off');
pxSize = 57.5;

%% simulation parameters
N = 48; % number of realizations per orientation
s = 400; % signal - median of fibril data
b = 15; % background - average of fibril data

thetaList = 10:10:80;
phiList = 0:15:165;
gammaT = 1; % rotational constraint - median of fibril data

img = nan(2*imgSz-1,2*imgSz-1,N);
loc_data_all = cell(length(thetaList),length(phiList));

for thetaInd = 1:length(thetaList)
    for phiInd = 1:length(phiList)
        thetaT = thetaList(thetaInd);
        phiT = phiList(phiInd);
        mux = sind(thetaT)*cosd(phiT);
        muy = -sind(thetaT)*sind(phiT);
        muz = cosd(thetaT);
        %% generate simulated images
        m = gammaT*[mux^2,muy^2,muz^2,mux*muy,mux*muz,muy*muz]+(1-gammaT)/3*[1,1,1,0,0,0];

        imgSim = zeros(2*imgSz-1);
        for i = 1:6
            imgSim = imgSim + B(:,:,i)*m(i);
        end

        for i = 1:N
            img(:,:,i) = poissrnd(imgSim*s+b);
        end

        %% RoSE-O estimation

        loc_data = [];
        obj.pixelUpsample = 1;
        obj.pixelSize = pxSize*2;

        % number of workers in parfor
        % parWorker = getenv('number_of_processors');
        parWorker = feature('numcores');
        for groupInd = 1:ceil(size(img,3)/parWorker)
            % RoSE-O works with images with odd size
            if groupInd ~= ceil(size(img,3)/parWorker) || rem(size(img,3),parWorker) == 0
                SMLM_img = img(:,:,groupInd*parWorker-parWorker+(1:parWorker));
            else
                SMLM_img = img(:,:,groupInd*parWorker-parWorker+(1:min([parWorker,rem(size(img,3),parWorker)])));
            end
            parfor subInd = 1:size(SMLM_img,3)
                [~,~,locDataTmp] = RoSEO_PolCam(obj, SMLM_img(:,:,subInd), b, FPSF, 'regVal', .25);
                if ~isempty(locDataTmp)
                    locDataTmp(:,1) = groupInd*parWorker+subInd-parWorker; % record frame number
                    for i = 1:size(locDataTmp,1)
                        % project to first order moments
                        [mux, muy, muz, rotMobil] = secondM2SymmConeWeighted(Bs,B_Ha,sumNorm,locDataTmp(i,5:10),locDataTmp(i,4),b);
                        locDataTmp(i,11:14) = [mux, -muy, muz, rotMobil]; % -muy to match coordinates
                    end
                    loc_data = [loc_data;locDataTmp];
                end
            end
        end
        disp(['[theta,phi] = [',num2str(thetaT),',',num2str(phiT),'] analyzed']);
        loc_data_all{thetaInd,phiInd} = loc_data;
    end
end

save('simResults.mat','loc_data_all');

%% visualizing estimation results

gamma = nan(N,length(thetaList),length(phiList));
theta = nan(N,length(thetaList),length(phiList));
phi = nan(N,length(thetaList),length(phiList));
locThreshold = pxSize*2;

for thetaInd = 1:length(thetaList)
    for phiInd = 1:length(phiList)
        loc = loc_data_all{thetaInd,phiInd};
        FN = 0; % false negative
        FP = 0; % false positive
        for i = 1:N
            locTmp = loc(loc(:,1)==i,:);
            if isempty(locTmp) || min(sqrt((locTmp(:,2)-pxSize*2).^2+(locTmp(:,3)-pxSize*2).^2)) > locThreshold % false negative
                FN = FN+1;
            else
                if size(locTmp,1) ~= 1 % false positive
                    [~,ind] = min(abs(locTmp(:,2)-pxSize*2)+abs(locTmp(:,2)-pxSize*2));
                    FP = FP+size(locTmp,1)-1;
                    locTmp = locTmp(ind,:);
                end
                % correct detection
                x(i,thetaInd,phiInd) = locTmp(2)-pxSize*2;
                y(i,thetaInd,phiInd) = locTmp(3)-pxSize*2;
                theta(i,thetaInd,phiInd) = acosd(abs(locTmp(13)));
                phiTmp = atan2d(locTmp(12),locTmp(11));
                if thetaInd == 7 && phiInd == 4
                    thetaInd;
                end
                if abs(phiTmp - phiList(phiInd)) >= 90
                    if phiTmp > phiList(phiInd)
                        phiTmp = phiTmp - 180;
                    else
                        phiTmp = phiTmp + 180;
                    end
                end
                phi(i,thetaInd,phiInd) = phiTmp;
                gamma(i,thetaInd,phiInd) = locTmp(14);
            end
        end
        false_positive(thetaInd,phiInd) = FP;
        false_negative(thetaInd,phiInd) = FN;
        
        % compute CRB
        [CRBr,CRBmu] = computeCRB(B,Bgradx,Bgrady,thetaList(thetaInd),phiList(phiInd),gammaT,s,b);
        thetaCRB(thetaInd,phiInd) = sqrt(CRBmu(1,1));
        phiCRB(thetaInd,phiInd) = sqrt(CRBmu(2,2));
        gammaCRB(thetaInd,phiInd) = sqrt(CRBmu(3,3));
        xCRB(thetaInd,phiInd) = sqrt(CRBr(1,1));
        yCRB(thetaInd,phiInd) = sqrt(CRBr(2,2));
    end
end

figure();
subplot(1,2,1);
imagesc(phiList,thetaList,false_positive/N*100);
colorbar; xlabel('\theta (deg.)'); ylabel('\phi (deg.)'); title('false positive %')
subplot(1,2,2);
imagesc(phiList,thetaList,false_negative/N*100);
colorbar; xlabel('\theta (deg.)'); ylabel('\phi (deg.)'); title('false negative %')
print(['simData_detection_signal=',num2str(s),'_background=',num2str(s),'_gamma=',num2str(gammaT)],'-dpng');

figure();
set(gcf,'position',[100,100,1800,800])
plotEst(1,squeeze(nanmean(x)),squeeze(nanstd(x)),xCRB,thetaList,phiList,{'x bias (nm)','\sigma_x (nm)','sqrt(CRB_x) (nm)'})
plotEst(2,squeeze(nanmean(y)),squeeze(nanstd(y)),yCRB,thetaList,phiList,{'y bias (nm)','\sigma_y (nm)','sqrt(CRB_y) (nm)'})
plotEst(3,squeeze(nanmean(theta))-thetaList',squeeze(nanstd(theta)),thetaCRB,thetaList,phiList,...
    {'\theta bias (deg.)','\sigma_{\theta} (deg.)','sqrt(CRB_{\theta}) (deg.)'})
plotEst(4,squeeze(nanmean(phi))-phiList,squeeze(nanstd(phi)),phiCRB,thetaList,phiList,...
    {'\phi bias (deg.)','\sigma_{\phi} (deg.)','sqrt(CRB_{\phi}) (deg.)'})
plotEst(5,squeeze(nanmean(gamma))-gammaT,squeeze(nanstd(gamma)),gammaCRB,thetaList,phiList,...
    {'\gamma bias','\sigma_{\gamma}','sqrt(CRB_{\gamma})'})
print(['simData_estimation_signal=',num2str(s),'_background=',num2str(s),'_gamma=',num2str(gammaT)],'-dpng');

function plotEst(i,estBias,estStd,estCRB,thetaList,phiList,figTitle)

subplot(3,5,i);
imagesc(phiList,thetaList,estBias); 
colorbar; xlabel('\theta (deg.)'); ylabel('\phi (deg.)'); title(figTitle{1}); 
set(gca,'xtick',0:45:180,'ytick',0:10:90);
subplot(3,5,i+5);
imagesc(phiList,thetaList,estStd); 
colorbar; xlabel('\theta (deg.)'); ylabel('\phi (deg.)'); caxis([min([estStd(:);estCRB(:)]),max([estStd(:);estCRB(:)])]); 
title(figTitle{2}); set(gca,'xtick',0:45:180,'ytick',0:10:90);
subplot(3,5,i+10);
imagesc(phiList,thetaList,estCRB); 
colorbar; xlabel('\theta (deg.)'); ylabel('\phi (deg.)'); caxis([min([estStd(:);estCRB(:)]),max([estStd(:);estCRB(:)])]);
title(figTitle{3}); set(gca,'xtick',0:45:180,'ytick',0:10:90);

end

