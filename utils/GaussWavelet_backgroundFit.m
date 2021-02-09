% 210208 revised by OZ

% 190604 - v3_1 TD - Changed saving format from uint16 to single. Also, use
% Tiff class for saving the floating point numbers.

% 190522 - v3 TD - Added "margin" in the fitting procedure

% 190226 - v2 TD - Analysis on multiple stacks. Some code tuning for
% compatibility with workstation (ll03) analysis

% Tianben Ding 190212
% This program estimate non-uniform background by fitting individual
% columns and raws by 1D 2 Gaussian. A wavelet filter is further applied to
% the sum.

function backgroundEst = GaussWavelet_backgroundFit(imgMean)

imgSz = size(imgMean,1);

sigmaMin = 40; % minimum sigma to be fit, pixel

% wavelet filtering para
sorh = 's';
wname = 'bior6.8';
level = 6;

colFit = nan(imgSz);
rowFit = nan(imgSz);

parfor i = 1:imgSz % scan from row#1 to row#end
    options = fitoptions('gauss2', 'Lower', [0 -Inf sigmaMin 0 -Inf sigmaMin]);
    f = fit((1:imgSz).',imgMean(i,:).','gauss2',options);
    a1 = f.a1;
    b1 = f.b1;
    c1 = f.c1;
    a2 = f.a2;
    b2 = f.b2;
    c2 = f.c2;
    rowFit(i,:) =  a1.*exp( -(((1:imgSz)-b1)./c1).^2 )  +  a2.*exp( -(((1:imgSz)-b2)./c2).^2 );
end
parfor i = 1:imgSz % scan from column#1 to column#end
    options = fitoptions('gauss2', 'Lower', [0 -Inf sigmaMin 0 -Inf sigmaMin]);
    f = fit((1:imgSz).',imgMean(:,i),'gauss2',options);
    a1 = f.a1;
    b1 = f.b1;
    c1 = f.c1;
    a2 = f.a2;
    b2 = f.b2;
    c2 = f.c2;
    colFit(:,i) =  a1.*exp( -(((1:imgSz)-b1)./c1).^2 )  +  a2.*exp( -(((1:imgSz)-b2)./c2).^2 );
end
ave1DFit = (rowFit + colFit)./2;

[C,S] = wavedec2(ave1DFit,level,wname);
thr = wthrmngr('dw2ddenoLVL','sqtwolog',C,S,'one');
backgroundEst = wdencmp('lvd',C,S,wname,level,thr,sorh);

end
