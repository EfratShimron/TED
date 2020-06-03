function  [dT_zoomed_Hybrid] = Hybrid(FullKspace,PARAMS,sampling_mask,t_jjj,lam)
% This function implements the 'K-space Hybrid Method' for temperature
% reconstruction, which was published in:
%     Gaur, Pooja, and William A. Grissom. "Accelerated MRI thermometry by direct
%     estimation of temperature from undersampled k?space data." MRM 73.5 (2015): 1914-1925.?
% This code is based on the Matlab toolbox provided by William Grissom in his
% website:
% https://vuiis.vumc.org/~grissowa/software.html



% ==========================================================
%                 K-space Hybrid Method Parameters
% ==========================================================
% From the help of kspace_hybrid_thermo.m:

%|  algp    Algorithm parameters structure containing (structure and each entry are optional):
%|              order       1             Polynomial order (default = 0)
%|              lam         [1 2]         l1 penalty weights for real and imaginary parts of m
%|                                        (default = 10^-6)
%|              beta        1             Roughness penalty weight for real
%|                                        and imaginary parts of m (default = 0)
%|                                        for second stage of algorithm. (radians; default = 0.01)
%|              dofigs      1             Display intermediate figures (default = 0)
%|              thiters     1             Number of CG iterations per theta update (default = 10)
%|              citers      1             Number of CG iterations per c update (default = 5)
%|              bls         1             (1) Use Boyd's backtracking line search (usually faster); (0) Use ours
%|              masknz      [Nx,Ny,Nz]    Mask of non-zero heating
%|                                        locations. This will cause
%|                                        the algorithm to skip the l1-regularized
%|                                        stage and go straight to the masked/unregularized stage
%|                                        (default = [])
%|

% The parameters of Gaur's method are with their DEFAULT values are:
% algp.maskthresh = -0.1/CONST; phase shift threshold for Step 6 (relaxation) (see their paper); equivalent to 0.1 deg C
% algp.order = 0; % polynomial order! was varied from 0 to 1 in version v9
% algp.lam = [10^-2 -1]; % sparsity regularization parameter
% algp.beta = 0; %2^-11; % roughness regularization parameter




disp('************ k-space hybrid method  *****************')


% ------ re-format data structure to fit the k-space method requirements ---------------
for t_ind = 1:max(PARAMS.t_rec_vec)
    for ncoil = 1:PARAMS.NC
        switch PARAMS.scanner_vendor
            case 'Philips'
                FullKspace_4kspace_method(:,:,t_ind,ncoil) = (squeeze(FullKspace(:,:,t_ind,ncoil)));

            otherwise
                FullKspace_4kspace_method(:,:,t_ind,ncoil) = conj(squeeze(FullKspace(:,:,t_ind,ncoil)));

                
        end
    end
end
data = permute(FullKspace_4kspace_method,[1 2 4 3]); % change order of dimensions to: Nx,Ny,NC,NT
if max(data(:))==0
    error('Kspace hybrid method - "data" array does not any contain data!')
    return
end

[Nx,Ny,Nc,Nt] = size(data); % # x,y-locs, coils, dynamics


% -------------- sampling mask --------------------
% to enable a fair comparison with other methods:
sampling_mask_4kspace_method = sampling_mask;
switch PARAMS.scanner_vendor 
    case 'Philips'
        sampling_mask_4kspace_method = imrotate(sampling_mask_4kspace_method,-90);
end

inds = find(sampling_mask_4kspace_method(:,1)==1);
if isempty(inds)==1
    error('K-space hybrid method: inds to sample is an empty vector')
    return
end

% ------ Create library of Baseline Images (for all coils) ----------
% recon the baselines
for ii = 1:Nc
   L(:,:,ii) = fftshift(ifft2(fftshift(squeeze(data(:,:,ii,PARAMS.t_baseline)))))*Nx*Ny; % added in v2, for the case of t_baseline~=1
end

% =========================== Reconstruction =========================
algp.lam = lam*[1 1];  % lam is in input of this function
algp.thiters = 10;
algp.citers = 5;
    
    
        
        t_ind = PARAMS.t_rec_vec(t_jjj);
                
        % sample horizontal lines (Efrat)
        dacc = permute(data(inds,:,:,t_ind),[3 1 2]); % extract sampled lines for this time frame
        dacc = dacc(:,:).';
        k = false(Nx,Nx);
        k(inds,:) = true;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % k-space recon from undersampled data with FFT's
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thetainit = zeros(Nx,Ny); % initialize temp phase shift map with zeros
        acqp.data = dacc; % accelerated data
        acqp.k = k(:,1); % k-space sampling mask - passing a vector tells ksh tp
        % premptively do FFT in fully-sampled dimension, and run in hybrid space
        
        
        acqp.L = L(:); % baseline 'library'
        algp.order = 0; % polynomial order        
        %algp.maskthresh = -0.05*PARAMS.CONST;
        algp.maskthresh = -0.1*PARAMS.CONST; % 
        
        
        algp.useGPU = false; % cast all variables to the GPU (Cartesian Only)
        algp.dofigs = 0; % show figures during the hybrid thermometry computation
        
        % ------- recon ----------
        [thetakacc,~,~,f,Ac] = kspace_hybrid_thermo(acqp,thetainit,algp);
        
        switch PARAMS.scanner_vendor
            case 'Philips'
                thetakacc_rotated_shifted = fftshift(thetakacc,2);
            otherwise
                thetakacc_rotated = imrotate(thetakacc,-180);
                thetakacc_rotated_shifted = circshift(thetakacc_rotated,1);
        end
        dT_Hybrid = real(thetakacc_rotated_shifted)/PARAMS.CONST; % the "real()" was taken from Grissoms's code
        dT_zoomed_Hybrid = dT_Hybrid(PARAMS.x1:PARAMS.x2,PARAMS.y1:PARAMS.y2);

