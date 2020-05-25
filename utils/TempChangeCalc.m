function [dT,dT_zoomed] = TempChangeCalc_v3(KspaceData_baseline,KspaceData_HIFU,S)
% xxx delete this note
% this function is based on the following function:
%[dT] = TempChangeCalc_1_frame_v2(KspaceData_baseline,KspaceData_HIFU,PARAMS)% newRealImages, newImgImages)

% This function computes dT, which is the Temperature Change
% related to the phase change between these two datasets:
% KspaceData_baseline, KspaceData_HIFU.
% First, we compute delta_Phi (the overall phase change) from the
% multi-coil K-space data, according to eq.[2] in the TED paper:
%                delta_Phi_t = atan( sum_i( (x_i_0)*conj(x_i_t)) ) where i=1:Nc is the coil index
% This calibrationless multi-coil merging method was introduced in:
% Bernstein et al., "Reconstructions of phase contrast, phased array multicoil data", MRM 1994;32(3)
% Second, we compute dT by dT = delta_Phi_t/PARAMS.CONST;

% KspaceData dimensions are NxNxNc

N = size(KspaceData_baseline,1);   % Each 2D slice is an NxN matrix
NC = size(KspaceData_baseline,3);  % number of coils

%-------------- Compute phase change ----------------
Imag_sum = zeros(N);
Real_sum = zeros(N);

for i=1:NC
    x_i_0 = ifft2c(squeeze(KspaceData_baseline(:,:,i)));        % coil i, baseline scan (time frame #1)
    x_i_t = ifft2c(squeeze(KspaceData_HIFU(:,:,i)));    % coil i, post-heating scan (in time frame t_ind)
    
    Re = real(x_i_0).*real(x_i_t) + imag(x_i_0).*imag(x_i_t);
    Im = imag(x_i_0).*real(x_i_t) - real(x_i_0).*imag(x_i_t);
    
    % This is a clinically-useful threshold:
    if S.PARAMS.combined_complex_th_flag~=0
        Re(find(Re<1))=1;
    end
    
    Imag_sum = Imag_sum + Im;
    Real_sum = Real_sum + Re;
end
delta_Phi_t = atan2(Imag_sum , Real_sum);       % this is an NxN matrix, obtained from the NC coils together

%-------------- compute temperature change -------------
dT = delta_Phi_t/S.PARAMS.CONST;  % delta_Phi_t is an NxN matrix

% --------- corrections needed for some vendors ----------
switch S.PARAMS.scanner_vendor
    case 'Philips'
        % here we make some simple corrections to the data due to the built-in
        % setup of Philips scanners, which invert the phase and perform fftshift
        % along a single dimension (not along both dimensions!).
        dT = -1*dT;
        dT = fftshift(dT,2);
end

dT_zoomed = dT(S.PARAMS.x1:S.PARAMS.x2,S.PARAMS.y1:S.PARAMS.y2);

end