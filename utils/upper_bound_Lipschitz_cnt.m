function l_max = upper_bound_Lipschitz_cnt(y, b, FPSF, recovStruct)
%upper_bound_Lipschitz_cnt computes an upper bound on the Lipschitz
%constant of the Poisson negative log-likelihood  TODO: cite the relevent paper
%->---
%input
%->---
%y
%N
%n_f
%b
%M
%up_sample
%---->-
%output
%---->-
%l_max

%x_channel
FXX = FPSF.FXX;
FYY = FPSF.FYY;
FZZ = FPSF.FZZ;
FXY = FPSF.FXY;
FXZ = FPSF.FXZ;
FYZ = FPSF.FYZ;

N = recovStruct.n_grid_p;
n_f = recovStruct.num_frames;
M = recovStruct.img_size;
upsample = recovStruct.upsample_factor;

% b_t=min(b);
l_max = (max(bsxfun(@times, y, 1 ./ b)).^2) .* max(reshape((fast_mul_fft(ones(6 * N, n_f))), M^2, n_f));


    function out_N1 = down_samp(x)
        out_N1 = x(1:upsample:end, 1:upsample:end, :);
    end
    function out_N2 = fast_mul_fft(x)

        out_N2 = [down_samp(real(ifft2(bsxfun(@times, FXX, fft2(xxgrid(x))) + bsxfun(@times, FYY, fft2(yygrid(x))) + ...
            bsxfun(@times, FZZ, fft2(zzgrid(x))) + bsxfun(@times, FXY, fft2(xygrid(x))) + bsxfun(@times, FXZ, fft2(xzgrid(x))) + ...
            bsxfun(@times, FYZ, fft2(yzgrid(x))))))];
       
        function out_N1_inN2 = xxgrid(x)
            out_N1_inN2_t = (reshape(x(1:N, :), sqrt(N), sqrt(N), n_f));
            out_N1_inN2 = (padarray(out_N1_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));

        end
        function out_N2_inN2 = yygrid(x)
            out_N2_inN2_t = (reshape(x(N + 1:2 * N, :), sqrt(N), sqrt(N), n_f));
            out_N2_inN2 = (padarray(out_N2_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end
        function out_N3_inN2 = zzgrid(x)
            out_N3_inN2_t = (reshape(x(2 * N + 1:3 * N, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end
        function out_N3_inN2 = xygrid(x)
            out_N3_inN2_t = (reshape(x(3 * N + 1:4 * N, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end
        function out_N3_inN2 = xzgrid(x)
            out_N3_inN2_t = (reshape(x(4 * N + 1:5 * N, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end

        function out_N3_inN2 = yzgrid(x)
            out_N3_inN2_t = (reshape(x(5 * N + 1:end, :), sqrt(N), sqrt(N), n_f));
            out_N3_inN2 = (padarray(out_N3_inN2_t(6:end - 5, 6:end - 5, :), [5, 5]));
        end

    end
end