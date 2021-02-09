function [FPSF,Bs,B_Ha,B,Bgradx,Bgrady,sumNorm] = computeBasis(PSFsz,PSFsz_small,vis)
%% prepare to generate basis images and gradient
% imaging parameters
NA = 1.4; n1 = 1.515; rm = NA/n1; lambda = 630e-9; pxSize = 57.5e-9; n2 = 1.33;

% sampling size
bfpRadius = 40;
bfpSampling = ceil(bfpRadius*2/rm);
if rem(bfpSampling,2) == 1
    bfpSampling = bfpSampling + 1;
end
N = ceil(lambda/pxSize/NA*bfpRadius);
if rem(N,2) == 1
    N = N + 1;
end

dx = n1*pxSize;
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),N),linspace(-1/(2*dx),1/(2*dx),N));
xBFP = lambda*eta;
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;

% mask representing lateral shift for gradient computation 
[X,Y] = meshgrid(1e-9*(1:N),1e-9*(1:N));
maskdx = exp(1j*2*pi*X/(N*pxSize));
maskdy = rot90(maskdx);
mask = cat(3,ones(N),maskdx,maskdy);

%% compute BFP electric field
theta1 = asin(rho);
theta2 = asin((n1/n2)*sin(theta1));
a = cos(theta2)./cos(theta1);
% Fresnel coefficients
ts = 2*n1./(n1+n2*a);
tpxy = 2*n1./(n1+n2./a);
tpz = 2*n1./(n2+n1*a);

Esx = -sin(phi).*ts;
Esy = cos(phi).*ts;
Epx = cos(phi).*cos(theta1).*tpxy;
Epy = sin(phi).*cos(theta1).*tpxy;
Epz = -sin(theta1).*(n1/n2).*tpz;

Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx); % Exx - Ex contributed by mux
Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy);
Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz);
Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx);
Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy);
Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz);

Exx(rho >= rho_max) = 0;
Exy(rho >= rho_max) = 0;
Exz(rho >= rho_max) = 0;
Eyx(rho >= rho_max) = 0;
Eyy(rho >= rho_max) = 0;
Eyz(rho >= rho_max) = 0;

Ex_bfp = cat(3,Exx,Exy,Exz);
Ey_bfp = cat(3,Eyx,Eyy,Eyz);

%% FT to find image plane electric field
if PSFsz*2 <= N
    PSFregion = N/2+(-PSFsz+1:PSFsz);
else
    PSFregion = PSFsz+(-N/2+1:N/2);
end

% detector polarizer grid
J = zeros(N,N,4);
J(1:2:end,1:2:end,[1,4]) = .5;
J(1:2:end,1:2:end,[2,3]) = -.5;
J(1:2:end,2:2:end,1) = 1;
J(2:2:end,1:2:end,4) = 1;
J(2:2:end,2:2:end,:) = .5;

for j = 1:3
    for i = 1:3
        Ex_bfpTmp(:,:,i) = Ex_bfp(:,:,i).*mask(:,:,j);
        Ey_bfpTmp(:,:,i) = Ey_bfp(:,:,i).*mask(:,:,j);
        ExTmp = fftshift(fft2(Ex_bfpTmp(:,:,i)));
        EyTmp = fftshift(fft2(Ey_bfpTmp(:,:,i)));
        Ex(:,:,i) = ExTmp.*J(:,:,1)+EyTmp.*J(:,:,2);
        Ey(:,:,i) = ExTmp.*J(:,:,3)+EyTmp.*J(:,:,4);
    end
    % match PSF size to image
    if PSFsz*2 <= N
        Elist{j} = [Ex(PSFregion,PSFregion,:),Ey(PSFregion,PSFregion,:)];
    else
        ExPad = zeros(PSFsz*2,PSFsz*2,3);
        EyPad = zeros(PSFsz*2,PSFsz*2,3);
        ExPad(PSFregion,PSFregion,:) = Ex;
        EyPad(PSFregion,PSFregion,:) = Ey;
        Elist{j} = [ExPad,EyPad];
    end
end

%% image plane basis
for j = 1:3
    E_tmp = Elist{j};
    for i = 1:3
        B_tmp(:,:,i) = abs(E_tmp(:,:,i)).^2;
    end
    B_tmp(:,:,4) = 2*real(E_tmp(:,:,1).*conj(E_tmp(:,:,2)));
    B_tmp(:,:,5) = 2*real(E_tmp(:,:,1).*conj(E_tmp(:,:,3)));
    B_tmp(:,:,6) = 2*real(E_tmp(:,:,2).*conj(E_tmp(:,:,3)));
    B_tmp = B_tmp/sum(sum(B_tmp(:,:,1)));
    Blist{j} = B_tmp(:,1:end/2,:)+B_tmp(:,end/2+1:end,:);
end
B = Blist{1};
Bgradx = Blist{2}-Blist{1};
Bgrady = Blist{3}-Blist{1};

B0 = B(PSFsz+[-PSFsz_small+1:PSFsz_small],PSFsz+[-PSFsz_small+1:PSFsz_small],:);
Bgradx0 = Bgradx(PSFsz+[-PSFsz_small+1:PSFsz_small],PSFsz+[-PSFsz_small+1:PSFsz_small],:);
Bgrady0 = Bgrady(PSFsz+[-PSFsz_small+1:PSFsz_small],PSFsz+[-PSFsz_small+1:PSFsz_small],:);

if rem(abs(PSFsz-N/2),2)
    B = [[B(2:2:end,2:2:end,:),B(2:2:end,1:2:end,:)];[B(1:2:end,2:2:end,:),B(1:2:end,1:2:end,:)]];
    Bgradx = [[Bgradx(2:2:end,2:2:end,:),Bgradx(2:2:end,1:2:end,:)];[Bgradx(1:2:end,2:2:end,:),Bgradx(1:2:end,1:2:end,:)]];
    Bgrady = [[Bgrady(2:2:end,2:2:end,:),Bgrady(2:2:end,1:2:end,:)];[Bgrady(1:2:end,2:2:end,:),Bgrady(1:2:end,1:2:end,:)]];
else
    B = [[B(1:2:end,1:2:end,:),B(1:2:end,2:2:end,:)];[B(2:2:end,1:2:end,:),B(2:2:end,2:2:end,:)]];
    Bgradx = [[Bgradx(1:2:end,1:2:end,:),Bgradx(1:2:end,2:2:end,:)];[Bgradx(2:2:end,1:2:end,:),Bgradx(2:2:end,2:2:end,:)]];
    Bgrady = [[Bgrady(1:2:end,1:2:end,:),Bgrady(1:2:end,2:2:end,:)];[Bgrady(2:2:end,1:2:end,:),Bgrady(2:2:end,2:2:end,:)]];
end

B = B(1:end-1,1:end-1,:); % RoSE-O is compatible with odd image size
Bgradx = Bgradx(1:end-1,1:end-1,:);
Bgrady = Bgrady(1:end-1,1:end-1,:);
    
% visualize basis images
if strcmp(vis,'on')
    figure(); set(gcf,'position',[100,100,1200,800])
    for i = 1:6
        subplot(2,3,i);
        imagesc(B(:,:,i)); axis image; set(gca,'xtick',[],'ytick',[]); colorbar;
    end
end
%% prepare basis structure for RoSE-O
FPSF.FXX = (fft2((fftshift(B(:,:,1)))));
FPSF.FYY = (fft2((fftshift(B(:,:,2)))));
FPSF.FZZ = (fft2((fftshift(B(:,:,3)))));
FPSF.FXY = (fft2((fftshift(B(:,:,4)))));
FPSF.FXZ = (fft2((fftshift(B(:,:,5)))));
FPSF.FYZ = (fft2((fftshift(B(:,:,6)))));

%gradients
FPSF.FXXdx = (fft2((fftshift(10^2*Bgradx(:,:,1)))));
FPSF.FXXdy = -(fft2((fftshift(10^2*Bgrady(:,:,1))))); % negative to match RoSE-O coordinates
FPSF.FYYdx = (fft2((fftshift(10^2*Bgradx(:,:,2)))));
FPSF.FYYdy = -(fft2((fftshift(10^2*Bgrady(:,:,2)))));
FPSF.FZZdx = (fft2((fftshift(10^2*Bgradx(:,:,3)))));
FPSF.FZZdy = -(fft2((fftshift(10^2*Bgrady(:,:,3)))));

% small PSFs for visualization and projection to first moment
Bs.XX = B0(:,:,1);
Bs.YY = B0(:,:,2);
Bs.ZZ = B0(:,:,3);
Bs.XY = B0(:,:,4);
Bs.XZ = B0(:,:,5);
Bs.YZ = B0(:,:,6);

sumNorm = sum(sum(Bs.XX));

% Hadamard products
B_Ha.aa = (Bs.XX) .* (Bs.XX);
B_Ha.ab = (Bs.XX) .* (Bs.YY);
B_Ha.ac = (Bs.XX) .* (Bs.ZZ);
B_Ha.ad = (Bs.XX) .* (Bs.XY);
B_Ha.ae = (Bs.XX) .* (Bs.XZ);
B_Ha.af = (Bs.XX) .* (Bs.YZ);

B_Ha.bb = (Bs.YY) .* (Bs.YY);
B_Ha.bc = (Bs.YY) .* (Bs.ZZ);
B_Ha.bd = (Bs.YY) .* (Bs.XY);
B_Ha.be = (Bs.YY) .* (Bs.XZ);
B_Ha.bf = (Bs.YY) .* (Bs.YZ);

B_Ha.cc = (Bs.ZZ) .* (Bs.ZZ);
B_Ha.cd = (Bs.ZZ) .* (Bs.XY);
B_Ha.ce = (Bs.ZZ) .* (Bs.XZ);
B_Ha.cf = (Bs.ZZ) .* (Bs.YZ);

B_Ha.dd = (Bs.XY) .* (Bs.XY);
B_Ha.de = (Bs.XY) .* (Bs.XZ);
B_Ha.df = (Bs.XY) .* (Bs.YZ);

B_Ha.ee = (Bs.XZ) .* (Bs.XZ);
B_Ha.ef = (Bs.XZ) .* (Bs.YZ);

B_Ha.ff = (Bs.YZ) .* (Bs.YZ);

end
