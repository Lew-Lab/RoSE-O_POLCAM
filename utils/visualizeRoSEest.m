function visualizeRoSEest(loc_data,sideL,binSz,figPrint)

locdata = loc_data(:,[2,3,11:14,4]);
locdata(locdata(:,5)<0,[3:5]) = -locdata(locdata(:,5)<0,[3:5]);

locdata(:,3) = atan2d(locdata(:,4),locdata(:,3));
locdata(:,4) = asind(locdata(:,5));
locdata(:,5) = locdata(:,6);
locdata(:,6) = locdata(:,7);
locdata(:,7) = [];

count = zeros(round(sideL*2/binSz+1));
phi = nan(round(sideL*2/binSz+1));
theta = nan(round(sideL*2/binSz+1));
gamma = nan(round(sideL*2/binSz+1));
brightness = nan(round(sideL*2/binSz+1));

indy = round((sideL+locdata(:,1))/binSz);
indx = round((sideL+locdata(:,2))/binSz);

for i = 1:sideL*2/binSz+1
    for j = 1:sideL*2/binSz+1
        if any((indx==i).*(indy==j))
            phiTmp = locdata(logical((indx==i).*(indy==j)),3);
            for k = 1:size(phiTmp)
                phiTmp(phiTmp-phiTmp(1)>90) = phiTmp(phiTmp-phiTmp(1)>90) - 180;
                phiTmp(phiTmp-phiTmp(1)<-90) = phiTmp(phiTmp-phiTmp(1)<-90) + 180;
            end
            thetaTmp = locdata(logical((indx==i).*(indy==j)),4);
            gammaTmp = locdata(logical((indx==i).*(indy==j)),5);
            brightnessTmp = locdata(logical((indx==i).*(indy==j)),6);
            phi(i,j) = nanmean(phiTmp.*brightnessTmp);
            theta(i,j) = nanmean(thetaTmp.*brightnessTmp);
            gamma(i,j) = nanmean(gammaTmp.*brightnessTmp);
            count(i,j) = length(phiTmp);
            brightness(i,j) = nanmean(brightnessTmp);
        end
    end
end

phi = phi./brightness;
theta = theta./brightness;
gamma = gamma./brightness;

phi(phi>180) = phi(phi>180)-180;
phi(phi<-180) = phi(phi<-180)+180;
phi(phi<0) = phi(phi<0) + 180;


phi(count<2) = nan;
theta(count<2) = nan;
gamma(count<2) = nan;

cMapP = parula(256);
cMapPhase = hsv(256);

for i = 1:sideL*2/binSz+1
    for j = 1:sideL*2/binSz+1
        if isnan(phi(i,j))
            phiPlot(i,j,:) = cat(3,0,0,0);
        else
            colorInd = round(phi(i,j)/180*(size(cMapPhase,1)-1))+1;
            phiPlot(i,j,:) = reshape(cMapPhase(colorInd,:),1,1,3);
        end
        if isnan(theta(i,j))
            thetaPlot(i,j,:) = cat(3,0,0,0);
        else
            colorInd = round(theta(i,j)/90*(size(cMapP,1)-1))+1;
            thetaPlot(i,j,:) = reshape(cMapP(colorInd,:),1,1,3);
        end
        if isnan(gamma(i,j))
            gammaPlot(i,j,:) = cat(3,0,0,0);
        else
            colorInd = round(gamma(i,j)*(size(cMapP,1)-1))+1;
            gammaPlot(i,j,:) = reshape(cMapP(colorInd,:),1,1,3);
        end
    end
end

figure(102);
set(gcf,'position',[100,100,1600,800]);
s1 = subplot(2,3,1);
countPlot = count;
countPlot(2:4,end-1-round(1000/binSz):end-2) = max(count(:));
imagesc(countPlot); set(gca,'ydir','normal');
colormap(s1,hot); set(gca,'xtick',[],'ytick',[]); colorbar; axis image; title('count'); caxis([0,prctile(count(count~=0),95)])
text(size(count,1)-3-round(1000/binSz),10,'1 μm','color','w','fontweight','bold')
subplot(2,3,2);
histogram(locdata(:,4));
xlabel('polar angle \theta (deg.)'); ylabel('count');
subplot(2,3,3);
histogram(locdata(:,5));
xlabel('\gamma (rotational constraint)'); ylabel('count');
s2 = subplot(2,3,4);
imagesc(phiPlot); set(gca,'ydir','normal');
colormap(s2,cMapPhase); caxis([0,180]); set(gca,'xtick',[],'ytick',[]); colorbar; axis image; title('azimuthal angle \phi (deg.)');
s3 = subplot(2,3,5);
imagesc(thetaPlot); set(gca,'ydir','normal');
colormap(s3,cMapP); caxis([0,90]); set(gca,'xtick',[],'ytick',[]); colorbar; axis image; title('polar angle \theta (deg.)');
s4 = subplot(2,3,6);
imagesc(gammaPlot); set(gca,'ydir','normal');
colormap(s4,cMapP); caxis([0,1]); set(gca,'xtick',[],'ytick',[]); colorbar; axis image; title('\gamma (rotational constraint)');
if strcmp(figPrint,'on')
    print('fibrils_NileRed','-dpng');
end

figure(103);
set(gcf,'position',[100,100,1000,900]);
plot((locdata(:,1)+100*[-cosd(locdata(:,3)),cosd(locdata(:,3))])'/1e3,(locdata(:,2)+100*[-sind(locdata(:,3)),sind(locdata(:,3))])'/1e3,'k');
axis image;
xlabel('x (μm)'); ylabel('y (μm)');
title('in-plane orientation');
if strcmp(figPrint,'on')
    print('fibrils_NileRed_phi','-dpng');
end

end