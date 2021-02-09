function [gammaf, recovStruct, loc_data] = RoSEO_PolCam(obj, SMLM_img, b, FPSF, varargin)
%RoSEO returns a list of molecular parameter estimates
%->---
%input
%->---
%SMLM_img :            array(m,m)    -stack of single-molecule...
%localization microscopy images(x-y channel concanotated)
%
%images are realizations of the following statistical measurement model:
%SMLM_img ~ Poisson(img+background),
%b:                    array(m,m) or array(1,1) -estimates of background for all
%frames
%---->-
%output
%---->-
%loc_data:                estimated molecular parameters per localization
%each row contains frame number, brightness, position(x,y) , second
%moments


%check for proper input images
%------------------------------------------------------------
if max(SMLM_img(:)) < 2 * mean(b(:))
    gammaf = [];
    loc_data = [];
    recovStruct = struct();
    return
end
if any(SMLM_img < 0)
    error('input image must not be background subtracted')
end

img_size = size(SMLM_img, 1); % image side-length size

if img_size ~= size(SMLM_img, 2)

    error('input image must be a square region')

end

%make sure image size is an odd number
if mode(img_size, 2) == 0
    error('side length (in number of camera pixels) of input image must be an odd value')
end

%Fourier transforms of the  point-spread function
%------------------------------------------------------------

%PSFs
FXX = FPSF.FXX;
FYY = FPSF.FYY;
FZZ = FPSF.FZZ;
FXY = FPSF.FXY;
FXZ = FPSF.FXZ;
FYZ = FPSF.FYZ;
%gradients
FXXdx = FPSF.FXXdx;
FXXdy = FPSF.FXXdy;
FYYdx = FPSF.FYYdx;
FYYdy = FPSF.FYYdy;
FZZdx = FPSF.FZZdx;
FZZdy = FPSF.FZZdy;

%break all the frames into a set of sub_frames of length subframe_l
%------------------------------------------------------------
s = opt2struct(varargin);

SMLM_img = reshape(SMLM_img, img_size^2, 1);

%validate b

if numel(b) == 1
    b = repmat(b, img_size^2, 1);
    b = reshape(b, img_size^2, 1);
elseif size(b, 1) == img_size || size(b, 2) == img_size
    b = reshape(b, img_size^2, 1);
else
    error('input background size does not match the input image.')
end

% fixed or learned parameters for recovery
%------------------------------------------------------------
upsample_factor = obj.pixelUpsample; %object space pixel upsampling
n_boundry_p = 2 + ceil(img_size/4); %number of camera pixels for removing boundry artifacts
pixel_size = obj.pixelSize;
r = 10^-2 * pixel_size / upsample_factor / 2;

n_grid_p = (upsample_factor * img_size - (upsample_factor - 1))^2; % number of grid points in the object space

% compute  lateral grid points
%----------------------------------------------------
pixel_size = obj.pixelSize;
lateral_grid_p = (-(img_size - 1) / 2:1 / upsample_factor:(img_size - 1) / 2) * pixel_size;

%smoothing parameter
%-------------------------------
mu = 1;

% joint sparsity regularizer (if brightness is higher than 5k
% photons, increase the reg_val

%if brightness is high use larger values e.g., 0.4
if isfield(s, 'regval')
    reg_val = s.regval;
else
    reg_val = .18;
end
% if single molecule level brightness reg_val=.13;

MaxIt1 = 200; %maximum number of iterations (first stage)

w = 1; %weights

%build recovery structure
%------------------------------------------------------------
recovStruct = struct();
recovStruct.img_size = img_size;
recovStruct.lateral_grid_p = lateral_grid_p;
recovStruct.num_frames = 1;
recovStruct.subframe_l = 1;
recovStruct.upsample_factor = upsample_factor;
recovStruct.n_grid_p = n_grid_p;
recovStruct.mu = mu;
recovStruct.reg_val = reg_val;
recovStruct.w = w;
%----------------------------------------------------
%molecular parameters consists of 6 joint brightness-orientation parameters and
%6 joint orientation-position parameter; total 12 parameter
gamma_init = zeros(12*n_grid_p, 1);

SMLM_img_r = zeros(sqrt(n_grid_p), sqrt(n_grid_p));

SMLM_img_r(1:upsample_factor:end, 1:upsample_factor:end, :) = ...
    reshape(SMLM_img(1:img_size^2), img_size, img_size); % up_sample

%gamma_init_t=A'*SMLM_img; computation in Fourier space.

gamma_init_t_xx = real(ifft2(bsxfun(@times, conj(FXX), fft2(SMLM_img_r))));
gamma_init_t_yy = real(ifft2(bsxfun(@times, conj(FYY), fft2(SMLM_img_r))));
gamma_init_t_zz = real(ifft2(bsxfun(@times, conj(FZZ), fft2(SMLM_img_r))));


gamma_init_t_xx = reshape(padarray(gamma_init_t_xx(n_boundry_p + 1:end - n_boundry_p, n_boundry_p + 1:end - n_boundry_p, :) ...
    , [n_boundry_p, n_boundry_p]), n_grid_p, 1);
gamma_init_t_yy = reshape(padarray(gamma_init_t_yy(n_boundry_p + 1:end - n_boundry_p, n_boundry_p + 1:end - n_boundry_p, :) ...
    , [n_boundry_p, n_boundry_p]), n_grid_p, 1);
gamma_init_t_zz = reshape(padarray(gamma_init_t_zz(n_boundry_p + 1:end - n_boundry_p, n_boundry_p + 1:end - n_boundry_p, :) ...
    , [n_boundry_p, n_boundry_p]), n_grid_p, 1);

gamma_init(1:n_grid_p, :) = repmat(sum(SMLM_img) ...
    ./sum(gamma_init_t_xx), n_grid_p, 1) .* (gamma_init_t_xx);

gamma_init(n_grid_p+1:2*n_grid_p, :) = repmat(sum(SMLM_img) ...
    ./sum(gamma_init_t_yy), n_grid_p, 1) .* (gamma_init_t_yy);

gamma_init(n_grid_p*2+1:3*n_grid_p, :) = repmat(sum(SMLM_img) ...
    ./sum(gamma_init_t_zz), n_grid_p, 1) .* (gamma_init_t_zz);


% upper bound on Lipschitz constant of the Poisson negative log-likelihood
%------------------------------------------------------------
l_max = upper_bound_Lipschitz_cnt(SMLM_img, b, FPSF, recovStruct);

l = l_max / 10;

% routiens
%------------------------------------------------------------
down_samp = @(x)x(1:upsample_factor:end, 1:upsample_factor:end, :);

xxgrid = @(x)(reshape(x(1:n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

yygrid = @(x)(reshape(x(n_grid_p + 1:2 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

zzgrid = @(x)(reshape(x(2 * n_grid_p + 1:3 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

xygrid = @(x)(reshape(x(3 * n_grid_p + 1:4 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

xzgrid = @(x)(reshape(x(4 * n_grid_p + 1:5 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

yzgrid = @(x)(reshape(x(5 * n_grid_p + 1:6 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

xxdxgrid = @(x)(reshape(x(6 * n_grid_p + 1:7 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

xxdygrid = @(x)(reshape(x(7 * n_grid_p + 1:8 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

yydxgrid = @(x)(reshape(x(8 * n_grid_p + 1:9 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

yydygrid = @(x)(reshape(x(9 * n_grid_p + 1:10 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

zzdxgrid = @(x)(reshape(x(10 * n_grid_p + 1:11 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

zzdygrid = @(x)(reshape(x(11 * n_grid_p + 1:12 * n_grid_p, :), sqrt(n_grid_p), sqrt(n_grid_p)));

%routines for computing Az via FFT

fA = @(x) abs(reshape(down_samp(real(ifft2(bsxfun(@times, FXX, fft2(xxgrid(x))) + ...
    bsxfun(@times, FYY, fft2(yygrid(x))) + bsxfun(@times, FZZ, fft2(zzgrid(x))) + ...
    bsxfun(@times, FXY, fft2(xygrid(x))) + bsxfun(@times, FXZ, fft2(xzgrid(x))) + ...
    bsxfun(@times, FYZ, fft2(yzgrid(x))) + bsxfun(@times, FXXdx, fft2(xxdxgrid(x))) + ...
    bsxfun(@times, FXXdy, fft2(xxdygrid(x))) + bsxfun(@times, FYYdx, fft2(yydxgrid(x))) + ...
    bsxfun(@times, FYYdy, fft2(yydygrid(x))) + bsxfun(@times, FZZdx, fft2(zzdxgrid(x))) + ...
    bsxfun(@times, FZZdy, fft2(zzdygrid(x)))))), img_size^2, 1));

% loop over sub_frames
%------------------------------------------------------------
%allocate space for first stage estimates

gamma_init = reshape(gamma_init, 12*n_grid_p, 1);

%prepare sub_frame structures

b_it = b;

SMLM_it = SMLM_img;
z = gamma_init;
gammaold = gamma_init;
w_it = w(:, :);
t = 1;
i = 1;
l_it = l(1);
recovStruct.w_it = w_it;
%sub_frame routines
%--------------------------------------------------------
%gradient of negative Poisson log-likelihood

gz = @(z)gradf2(z, SMLM_it, b_it, FPSF, recovStruct);

zq = @(z, gz_var, l_var)projection_cone2(z-repmat((1./l_var), 12 * n_grid_p, 1).*gz_var, r, recovStruct);

%                 zq=@(z,gz_var,l_var)proxmu(z-repmat((1./l_var),6*n_grid_p,1).*gz_var,recovStruct);


%Poisson negative log-likelihood

f1 = @(var)sum(var) - sum(SMLM_it.*log(bsxfun(@plus, var, b_it)));

q1 = @(var, z, gz_var, zq_var, l_var)sum(var) - sum(SMLM_it.*log(bsxfun(@plus, var, b_it))) + sum((gz_var).*(zq_var - z), 1) + ...
    (l_var / 2) .* sum((z-zq_var).^2, 1);

%handle for updating gamma via proximal operator

gammanew_update = @(z, gz_var, l_var) projection_cone2(z-repmat((1./l_var), 12 * n_grid_p, 1).* ...
    (gz_var + gradhmu2(z, recovStruct)), r, recovStruct);

%                 gammanew_update=@(z,gz_var,l_var) proxmu( z-repmat((1./l_var),6*n_grid_p,1).*...
%                     gz_var,recovStruct);


% set boundry condition - added by OZ for 4-ch pol camera
z0 = ones(img_size,img_size,12);
z0([1:n_boundry_p,end-n_boundry_p:end],:,:) = 0;
z0(:,[1:n_boundry_p,end-n_boundry_p:end],:) = 0;
z0 = reshape(z0,[],1);

while (i < MaxIt1)

    if (mod(i, 20) == 0 && i < 250)
        l_it = l_it / 5;
    end

    %backtracking line search
    k = 1;
    eta = 1.1; % line-search parameter
    gz_it = gz(z);
    zq_it = zq(z, gz_it, l_it);
    Azq = fA((zq_it));
    Az = fA((z));
    comp1 = f1(Azq) > (q1(Az, z, gz_it, zq_it, l_it) + .05); % descent condition

    while (any(comp1))

        l_it(comp1) = (eta(comp1).^k) .* l_it(comp1);

        zq_it = zq(z, gz_it, l_it);

        Azq = fA((zq_it));

        comp1 = f1(Azq) > (q1(Az, z, gz_it, zq_it, l_it) + .05);
        %disp(f1(Azq)-(q1(Az,z,gz_it,zq_it,l_it)))
        %disp(l_it)
        k = k + 1;
    end

    %update estimates
    %-----------------------------------------------------

    ll = l_it + 1 / mu; %Lipschitz constant of the smoothed objectvie function

    z = z.*z0;
    
    gammanew = gammanew_update(z, gz_it, ll);

    tnew = 1 / 2 + (sqrt(1 + 4 * t^2)) / 2;

    %momentum term of FISTA

    z = gammanew + ((t - 1) / tnew) * (gammanew - gammaold);

    gammaold = gammanew; % record the most recent estimate

    t = tnew;
    i = i + 1;

end

%record final estimates
gammahat1 = gammanew;
%num_char=progress_bar(round(nn/(n_sub_frame+1),2),num_char,20);


%===========================

%% GradMap copmutation

%===========================
[gradMap, gammhat2] = grad_map2(gammahat1, recovStruct);
if ~any(gammhat2)
    gammaf = [];
    loc_data = [];
    return
end
%===========================

%% adaptive, constrained  maximum likelihood

%===========================

%-----------------------------------------------------

%% summary

% After identifying  number of emitters correctly,
%we proceede to accurately and precisely estimate
%emitter's position and second moments
%We assume each channel has distinct offset along x and y
%step 1- transform the initial estimate to the
%new representation
%step2- design layers of the optimization program
%step 3- design a mapping from joint space to the
%molecular space
%-----------------------------------------------------

%% step 1- mapping initial estimates to

%representation with distinct channel offsets

% SKIPPING FOR NOW

%% step 2-

%fixed recovery parameters

MaxIt2 = 120;
eta = 1.2; % backtracking exponentiation parameter

%------------------------------------------------------------

% the complement of support set
indx_xx = (gammhat2(1:n_grid_p) == 0);
indx_yy = (gammhat2(n_grid_p + 1:2 * n_grid_p) == 0);
indx_zz = (gammhat2(2 * n_grid_p + 1:3 * n_grid_p) == 0);

indx = (indx_xx .* indx_yy .* indx_zz) > 0;


%
b_it = b;
z = gammhat2(:, :);
gamma_old = gammhat2(:, :);
SMLM_it = SMLM_img(:, :);
l_it = l(1) / 100000;

%set iteration parameters
t = 1;
i = 1;

%sub_frame routines
%---------------------------------------------------------

%gradient of negative Poisson log-likelihood

gz = @(z)gradf2(z, SMLM_it, b_it, FPSF, recovStruct);

zq = @(z, gz_var, l_var, indx_xx, indx_yy, indx_zz)projection_cone2(z-repmat((1./l_var), 12 * n_grid_p, 1).*gz_var, ...
    3*r, recovStruct, indx_xx, indx_yy, indx_zz);

%Poisson negative log-likelihood

f1 = @(var)sum(var) - sum(SMLM_it.*log(bsxfun(@plus, var, b_it)));

q1 = @(var, z, gz_var, zq_var, l_var)sum(var) - sum(SMLM_it.*log(bsxfun(@plus, var, b_it))) + sum((gz_var).*(zq_var - z), 1) + ...
    (l_var / 2) .* sum((z-zq_var).^2, 1);

%handle for updating gamma via proximal operator

gammanew_update = @(z, gz_var, l_var, indx_xx, indx_yy, indx_zz) projection_cone2(z-repmat((1./l_var), 12 * n_grid_p, 1).* ...
    (gz_var), 3*r, recovStruct, indx_xx, indx_yy, indx_zz);

while (i < MaxIt2)

    k = 1;
    gz_it = gz(z);
    zq_it = zq(z, gz_it, l_it, indx, indx, indx);
    Azq = fA(zq_it);
    Az = fA(z);
    comp1 = f1(Azq) > (q1(Az, z, gz_it, zq_it, l_it) + .05);
    %backtracking line search

    while (any(comp1))

        l_it(comp1) = (eta^k) * l_it(comp1);
        zq_it = zq(z, gz_it, l_it, indx, indx, indx);
        Azq = fA(zq_it);
        comp1 = f1(Azq) > (q1(Az, z, gz_it, zq_it, l_it) + .05);
        k = k + 1;
    end

    %update estimates
    %-----------------------------------------------------
    gammanew = gammanew_update(z, gz_it, l_it, indx, indx, indx);
    tnew = 1 / 2 + (sqrt(1 + 4 * t^2)) / 2;
    %momentum term of FISTA
    z = gammanew + ((t - 1) / tnew) * (gammanew - gamma_old);
    %record current estimate
    gamma_old = gammanew;
    t = tnew;
    i = i + 1;
end
%record final estimates
gammaf = gammanew;


%% mapping the joint parameters back to

% second moments and position

loc_data = get_loc_data2(gammaf, recovStruct);
