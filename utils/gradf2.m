function g = gradf2(gamma, SMLM_img, b, FPSF, recovStruct)
%gradf computes the gradient of the negative Poisson log-likelihood
%w.r.t. molecular parameters at the current estimate gamma
%->-----
% input
%->-----
% gamma:        array (N*3,n_f) -current estimate of the molecular parameters
% SMLM_img:     array (M,n_f)   -input camera images
% b:            array (M,n_f)   -background

% N:                            -number of grid points
% M:                            -image size in number of camera pixels
% n_f:                          -number of frames
% upsamplef:                    -object space upsample factor
%----->-
% output
%----->-
% g:                            -array (N*3,n_f) gradient of the negative Poisson
%                                log-likelihood w.r.t. gamma

%% 1- global parameters

N = recovStruct.n_grid_p;
M = recovStruct.img_size;
n_f = recovStruct.subframe_l;
upsamplef = recovStruct.upsample_factor;

%number of boundry pixels to guard against artifacts
n_boundry_p = 5;

%% 2- extract x-y channel parameters

%------------------------------------------------------------
%x_channel Fourier transfor of PSFs
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

boundxryIndx = [1:n_boundry_p, sqrt(N) - n_boundry_p + 1:sqrt(N)];

%% 3- computations are performed in Fourier domain

%------------------------------------------------------------
c_temp = (down_sample(fast_mul_fft(reshapeMat(gamma))));

c = bsxfun(@plus, reshape(c_temp, M^2, n_f), b);

g = fast_transp_mul_fft_x((1-SMLM_img./c));

%% utility functions

%------------------------------------------------------------


    function out = reshapeMat(A)
        out = reshape(A, sqrt(N), sqrt(N), 12, n_f);
    end

    function down_sampled = down_sample(x)
        down_sampled = x(1:upsamplef:end, 1:upsamplef:end, :);
    end


    function out_N = fast_mul_fft(x)

        out_N = real(ifft2(bsxfun(@times, FXX, fft2(xxgrid(x))) + bsxfun(@times, FYY, fft2(yygrid(x))) + ...
            bsxfun(@times, FZZ, fft2(zzgrid(x))) + bsxfun(@times, FXY, fft2(xygrid(x))) + bsxfun(@times, FXZ, fft2(xzgrid(x))) + ...
            bsxfun(@times, FYZ, fft2(yzgrid(x))) + bsxfun(@times, FXXdx, fft2(xxdxgrid(x))) + bsxfun(@times, FXXdy, fft2(xxdygrid(x))) + ...
            bsxfun(@times, FYYdx, fft2(yydxgrid(x))) + bsxfun(@times, FYYdy, fft2(yydygrid(x))) + bsxfun(@times, FZZdx, fft2(zzdxgrid(x))) + ...
            bsxfun(@times, FZZdy, fft2(zzdygrid(x)))));
    end


    function out_N1_inN2_t = xxgrid(x)
        out_N1_inN2_t = (reshape(x(:, :, 1, :), sqrt(N), sqrt(N), n_f));
        %         out_N1_inN2=(padarray(out_N1_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N1_inN2_t(boundxryIndx, :) = 0;
        out_N1_inN2_t(:, boundxryIndx) = 0;

    end

    function out_N2_inN2_t = yygrid(x)
        out_N2_inN2_t = (reshape(x(:, :, 2, :), sqrt(N), sqrt(N), n_f));
        %         out_N2_inN2=(padarray(out_N2_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N2_inN2_t(boundxryIndx, :) = 0;
        out_N2_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = zzgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 3, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = xygrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 4, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = xzgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 5, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = yzgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 6, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = xxdxgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 7, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2_t=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end
    function out_N3_inN2_t = xxdygrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 8, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = yydxgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 9, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2_t=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = yydygrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 10, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3_inN2_t = zzdxgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 11, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end
    function out_N3_inN2_t = zzdygrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 12, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3 = fast_transp_mul_fft_x(z)
        z_up = up_samp(z);
        fz = fft2(z_up);
        xx_temp = ifft2(bsxfun(@times, conj(FXX), fz));
        yy_temp = ifft2(bsxfun(@times, conj(FYY), fz));
        zz_temp = ifft2(bsxfun(@times, conj(FZZ), fz));
        xy_temp = ifft2(bsxfun(@times, conj(FXY), fz));
        xz_temp = ifft2(bsxfun(@times, conj(FXZ), fz));
        yz_temp = ifft2(bsxfun(@times, conj(FYZ), fz));
        xxdx_temp = ifft2(bsxfun(@times, conj(FXXdx), fz));
        xxdy_temp = ifft2(bsxfun(@times, conj(FXXdy), fz));
        yydx_temp = ifft2(bsxfun(@times, conj(FYYdx), fz));
        yydy_temp = ifft2(bsxfun(@times, conj(FYYdy), fz));
        zzdx_temp = ifft2(bsxfun(@times, conj(FZZdx), fz));
        zzdy_temp = ifft2(bsxfun(@times, conj(FZZdy), fz));

        indx = [1:n_boundry_p, size(xx_temp, 1) - n_boundry_p + 1:size(xx_temp, 1)];
        xx_temp(indx, :) = 0;
        xx_temp(:, indx) = 0;
        yy_temp(indx, :) = 0;
        yy_temp(:, indx) = 0;
        zz_temp(indx, :) = 0;
        zz_temp(:, indx) = 0;
        xy_temp(indx, :) = 0;
        xy_temp(:, indx) = 0;
        xz_temp(indx, :) = 0;
        xz_temp(:, indx) = 0;
        yz_temp(indx, :) = 0;
        yz_temp(:, indx) = 0;
        xxdx_temp(:, indx) = 0;
        xxdy_temp(:, indx) = 0;
        yydx_temp(:, indx) = 0;
        yydy_temp(:, indx) = 0;
        zzdx_temp(:, indx) = 0;
        zzdy_temp(:, indx) = 0;


        out_N3 = real([reshape(xx_temp, N, n_f); ...
            reshape(yy_temp, N, n_f); ...
            reshape(zz_temp, N, n_f); ...
            reshape(xy_temp, N, n_f); ...
            reshape(xz_temp, N, n_f); ...
            reshape(yz_temp, N, n_f); ...
            reshape(xxdx_temp, N, n_f); ...
            reshape(xxdy_temp, N, n_f); ...
            reshape(yydx_temp, N, n_f); ...
            reshape(yydy_temp, N, n_f); ...
            reshape(zzdx_temp, N, n_f); ...
            reshape(zzdy_temp, N, n_f); ...
            ]);


        function out_N1_N3 = up_samp(x)
            out_N1_N3 = zeros(sqrt(N), sqrt(N), n_f);
            out_N1_N3(1:upsamplef:end, 1:upsamplef:end, :) = reshape(x, M, M, n_f);
        end


    end

end