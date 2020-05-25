function [dT] = TempChangeCalc_1_frame_v2(KspaceData_baseline,KspaceData_HIFU,PARAMS)% newRealImages, newImgImages)
% This function first computes delta_Phi (the overall phase change) from the
% multi-coil K-space data, according to eq.[2] in the TED paper:
%                delta_Phi_t = atan( sum_i( (x_i_0)*conj(x_i_t)) ) where i=1:Nc is the coil index             
% This calibrationless multi-coil merging method was introduced in:
% Bernstein et al., "Reconstructions of phase contrast, phased array multicoil data", MRM 1994;32(3)


% KspaceData dimensions are: NxNxNc

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
    if PARAMS.combined_complex_th_flag~=0
        Re(find(Re<1))=1;
    end
        
    Imag_sum = Imag_sum + Im;
    Real_sum = Real_sum + Re;
end
delta_Phi_t = atan2(Imag_sum , Real_sum);       % this is an NxN matrix, obtained from the NC coils together

%-------------- compute temperature change -------------
dT = delta_Phi_t/PARAMS.CONST;  % delta_Phi_t is an NxN matrix



