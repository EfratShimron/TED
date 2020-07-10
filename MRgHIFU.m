classdef MRgHIFU
       
    properties
        demo
        PARAMS
        FLAGS
        
        FullKspace_all_slices
        FullKspace  % fully sampled k-space data of a specific slice
        
        sampling_mask
        sampling_mask_with_ACS
        
        GOP
        
        dT_gold
        dT_gold_zoomed
        
    end
    
    
    
    
    methods
        
        %% =========== Initialization ==============
        function S = MRgHIFU(demo,R_vec) % Initialization
            
            S.demo = demo;
            S=set_demo_params(S);
            S.PARAMS.R_vec = R_vec;
            
        end
        
        %% =========== set params ===================
        function S=set_demo_params(S)
            % =============== Initialize demo-specific Parameters ===================
            
            PARAMS.t_baseline = 1;  % for all demos, except the ones in which a different baseline frame is defined
            
            switch S.demo
                
                case 'Gel_phantom_demo'  % Insightec gel phantom, series 25
                    PARAMS.title = 'Gel phantom GE';
                    PARAMS.scanner_vendor = 'GE';
                    
                    % Acquisition params
                    PARAMS.B0 = 3; % [Tesla].
                    PARAMS.TE = 11.2;
                    PARAMS.N = 256; % Image size is NxN

                    % Reconstruction params
                    PARAMS.t_example = 5; % a time frame for plotting examples (for fully-sampled data)
                    PARAMS.t_rec_vec = 5; % time frames for dT reconstruction
              
                    PARAMS.wavWeight_SPIRIT = 0.00010;  % typical values for grid search:  %[1e-7 1e-6 1e-5 1e-4 0.0005 1e-3 0.0025 0.005 1e-2 0.025 0.05 ]; 
                    PARAMS.wavWeight_TED = 0.00045; % % typical values for grid search:  %[1e-7 1e-6 1e-5 1e-4 0.0005 1e-3 0.0025 0.005 1e-2 0.025 0.05 ]; 
                    PARAMS.nIter = 30; % number of Compressed Sensing iterations in TED and l1-SPIRiT
                    PARAMS.lam = 0.33*1e-3; % for k-space hybrid method
                    
                    PARAMS.combined_complex_th_flag = 0;  % This flag is for the multi-coil phase merging. It should be 1 only for in-vivo data
                   
                    % location of k-space data
                    PARAMS.kspace_data_filename = 'MRgHIFU_data/Gel_phantom_kspace_data'
                    
                   % plotting params
                    PARAMS.cmin = -25; PARAMS.cmax = 25;    % color limits for temperature plots (in degree Celsius)
                    PARAMS.x1=95; PARAMS.x2=155; PARAMS.y1=120;  PARAMS.y2=146;  % zoom-in coordinates - tight area around the HIFU heated zone.
                    PARAMS.mag_crop_x = 16:200;
                    PARAMS.mag_crop_y = 60:185;
                    PARAMS.mag_max = 0.04;
                    PARAMS.slices_for_recon = 1;                   
                    
                case 'Agar_phantom_demo' % Grissom's data
                    PARAMS.title = 'Agar phantom';
                    PARAMS.scanner_vendor = 'Philips';
                    
                    % Acquisition params
                    PARAMS.B0 = 3; % [Tesla].
                    PARAMS.TE = 12; %
                    PARAMS.N = 128; % Image size is NxN
                    
                    % Reconstruction params
                    PARAMS.t_example = 20;         % a time frame for plotting examples (for fully-sampled data)
                    PARAMS.t_rec_vec = [18 20 22]; % time frames for dT reconstruction  (other time frames can be chosen)

                    PARAMS.wavWeight_SPIRIT = 0.01;%  % typical values for grid search:  %[1e-7 1e-6 1e-5 1e-4 0.0005 1e-3 0.0025 0.005 1e-2 0.025 0.05 ]; 
                    PARAMS.wavWeight_TED = 0.005; %% typical values for grid search:  %[1e-7 1e-6 1e-5 1e-4 0.0005 1e-3 0.0025 0.005 1e-2 0.025 0.05 ]; 
                    PARAMS.nIter = 30; %  number of Compressed Sensing iterations in TED and l1-SPIRiT
                    PARAMS.lam = 1e-6; % for k-space hybrid method                                                    
                    
                    PARAMS.combined_complex_th_flag = 0;  % This flag is for the multi-coil phase merging. It should be 1 only for in-vivo data                    
                    
                    % plotting params
                    PARAMS.cmin = -30; PARAMS.cmax = 30;    % color limits for temperature plots (in degree Celsius)
                    PARAMS.mag_max = 5;                     % upper color limit for magnitude plots
                    PARAMS.x1=55;  PARAMS.x2=75;  PARAMS.y1=45;  PARAMS.y2=65; % zoom-in coordinates
                    
                    % location of k-space data
                    PARAMS.kspace_data_filename = 'MRgHIFU_data/Agar_phantom_kspace_data';
  
                    % plotting params
                    PARAMS.mag_crop_x = 40:90;
                    PARAMS.mag_crop_y = 10:105;
                    PARAMS.mag_max = 0.06;
                    PARAMS.slices_for_recon = 1;
                    
                    
            end
            
            % =========================== General Parameters =========================
            PARAMS.alpha = 9.0900e-06;
            PARAMS.gamma = 42.58;
            PARAMS.CONST = -2*pi*PARAMS.B0*PARAMS.gamma*PARAMS.alpha*PARAMS.TE ;%
            
            % ==================== HIFU mask ===================
            % The HIFU mask has 1s inside a block around the heated area, and 0
            % outside. It should be defined manually (or set to 1 everywhere)
            PARAMS.HIFU_MASK = zeros(PARAMS.N);
            PARAMS.HIFU_MASK(PARAMS.x1:PARAMS.x2,PARAMS.y1:PARAMS.y2) = 1;
            
            % ============== calibration params ==============
            
            switch S.demo
                case 'Agar_phantom_demo'
                    PARAMS.kSize = [3,3] %[3,3];  % SPIRiT kernel size
                    PARAMS.CalibTyk = 0.01;
                otherwise
                    PARAMS.kSize = [5,5];  % SPIRiT kernel size
                    PARAMS.CalibTyk = 0.01;  % Tykhonov regularization in the calibration
                    
            end
            
            S.PARAMS = PARAMS;
            
        end
        
        % ====================================================================================
        
        function S = load_kspace_data(S)
            disp('loading fully sampled data')
            
            load(S.PARAMS.kspace_data_filename)

            % The matrix FullKspace_all_slices has the following
            % dimensions: [N, N, NT, NC, NS]
            % where
            % N = image size
            % NT = number of Time-frames that were acquired during MRgHIFU session
            % NC = Number of Coils
            % NS = Number of Slices
            S.FullKspace_all_slices = FullKspace;            
            
            
            [NX,NY,NT,NC,NS]= size(S.FullKspace_all_slices);
            S.PARAMS.N = NX;  % assuming that NX=NY, we define them as N. If you have a non-square matrix, feel free to modify this.
            S.PARAMS.NT = NT;
            S.PARAMS.NC = NC;
            S.PARAMS.NS = NS;
            
            clear FullKspace
            
        end
        
        
        % ================= load sampling mask ======================
        function S = load_sampling_mask(S,R)
            % ======================= load/create sampling mask ============================
            load(['VarDensSampling/Sampling_mask_1D_var_dens_R',num2str(R),'_N',num2str(S.PARAMS.N),'.mat'])
            
            switch S.PARAMS.scanner_vendor % added in version 5
                case 'GE'
                    sampling_mask = imrotate(sampling_mask,90);  % horizontal lines (not columns) should be sampled for compatibility with the k-space method
            end
            
            % add a fully sampled center to the mask (needed for SPIRIT/GRAPPA kernel calibration)
            N = size(sampling_mask,1);
            mask_with_ACS = sampling_mask; % ACS = AutoCalibration Signal
            N_mid = round(N/2)+1;
             if N==256             
            vec = (N_mid-11):1:(N_mid+11);
            elseif N==128
                vec = (N_mid-15):1:(N_mid+15);
            end
            mask_with_ACS(vec,vec)=1;
            
            S.sampling_mask = sampling_mask;
            S.sampling_mask_with_ACS = mask_with_ACS;
            
            [CalibSize, dcomp] = getCalibSize(S.sampling_mask_with_ACS);
            S.PARAMS.dcomp = dcomp; % density compesnation matrix. see more explanation inside the funciton getCalibsize
            S.PARAMS.CalibSize = CalibSize; % see more explanation inside the funciton getCalibsize
        end
        
        % =====================================================
        function S = scale_and_calc_GOP(S)
            % This function computes the following things for both TED and
            % l1-SPIRiT:
            % - compensation mask (DATAcomp)
            % - scaline factor (scale_fctr) 
            % - the SPIRiT kernels (GOP)
            
            % =================== scale k-space data =================
            % Notes: 
            % 1. scaling the data such that the zero-filled density-compensated     
            % k-space norm is 1 is useful in order to use similar                       
            % regularization penalty values for different problems.                         
            %
            % 2. REMEMBER - for MR-guided-HIFU, temperature change is linearly
            % dependent upon the phase change, hence the data must be re-scaled back to
            % its original scale before the temperature change is computed.            
            
            % 3. to compute DATAcom, scale_fctr and the SPIRiT kernels, we sample the baseline data (i.e. data of time frame #1), 
            % since these factors and kernels will be used for reconstruction in later time frames. 

            DATA_baseline = squeeze(S.FullKspace(:,:,S.PARAMS.t_baseline,:));
            
            DATA_baseline_sampled = DATA_baseline.*repmat(S.sampling_mask_with_ACS,[1,1,S.PARAMS.NC]); % multiply with sampling matrix            
            
            % compute compensation mask & scaling factor
            DATAcomp = DATA_baseline_sampled.*repmat(S.PARAMS.dcomp,[1,1,S.PARAMS.NC]);  % sampled data x density compensation matrix
            scale_fctr = norm(DATAcomp(:))/sqrt(S.PARAMS.NC)/1;
            
            DATA_baseline_sampled = DATA_baseline_sampled/scale_fctr;
           
            % Calibrate SPIRiT kernels (for SPIRiT & TED)
            kCalib = crop(DATA_baseline_sampled,[S.PARAMS.CalibSize,S.PARAMS.NC]);
            kernel = zeros([S.PARAMS.kSize,S.PARAMS.NC,S.PARAMS.NC]);
            [AtA] = corrMatrix(kCalib,S.PARAMS.kSize);
            for n=1:S.PARAMS.NC
                kernel(:,:,:,n) = calibrate(AtA,S.PARAMS.kSize,S.PARAMS.NC,n,S.PARAMS.CalibTyk);
            end
            GOP = SPIRiT(kernel, 'fft',[S.PARAMS.N,S.PARAMS.N]); % Create the SPIRiT calibration kernels 
            
            S.PARAMS.scale_fctr = scale_fctr;
            S.GOP = GOP;
            
        end
        
        
        % ===================================================================
        function [dT_rec,dT_rec_zoomed,dT_NRMSE_per_iter,dT_RMSE_per_iter,mag_NRMSE_per_iter,Xsqr_gold,Xsqr_rec] = TED_and_SPIRiT(data, S, wavWeight,show,DATA_full_samp,DATA_baseline,R,calc_error_per_iter_flag,method_flag)
            % This function computes temperature reconstruction using two methods:
            % 1. l1-SPIRIT
            % 2. Temporal Differences (TED)
            
            % -----------------------------------------------------------------------
            %                       Documentation from ESPIRIT toolbox
            %
            % res = pocsSPIRIT(y, GOP, nIter, x0, wavWeight, show)
            %
            % Implementation of the Cartesian, POCS l1-SPIRiT reconstruction
            %
            % Input:
            %		y = data = Undersampled k-space data. Make sure that empty entries are zero
            %			       or else they will not be filled. (This is the zero-filled scaled k-space data)
            %		GOP -	the SPIRiT kernel operator obtained by calibration
            %		nIter -	Maximum number of iterations
            %		x0 - initial guess (in k-space, i.e. this is the zero-filled scaled k-space data)
            %		wavWeight - wavlet threshodling parameter
            % 		show - >1 to display progress (slower)
            %
            % Outputs:
            %		res - Full k-space matrix
            %
            % (c) Michael Lustig 2007
            %
            
            
            
            mask = (data==0);
            
            % ----------- calc gold standard magnitude image ----------------
            X_gold_baseline = ifft2c(DATA_baseline);
            
            X_gold = ifft2c(DATA_full_samp);
            Xsqr_gold = sqrt(sum(abs(X_gold).^2,3));  % for later computations of temperature RMSE
            [Xsqr_gold] = mag_corrections_per_vendor(Xsqr_gold,S.PARAMS);
            
          
            % ==========================================================
            
            % find the closest diadic size for the images
            [sx,sy,nc] = size(data);
            ssx = 2^ceil(log2(sx));
            ssy = 2^ceil(log2(sy));
            ss = max(ssx, ssy);
            W = Wavelet('Daubechies',4,4);
            
            x0 = data;  % initial guess in k-space = sampled k-space data (zero-filled)
            
            % % ------ show initial guess  ------
            % if show==1
            %     X = ifft2c(x0);
            %     Xsqr = sqrt(sum(abs(X).^2,3));
            %     figure; imshow(Xsqr,[],'InitialMagnification',400);
            %     title(['initial guess']) 
            %     drawnow
            % end
            
            switch method_flag
                case 'SPIRIT'
                    wavWeight_SPIRIT = wavWeight;
                    % =====================================================================
                    %                           l1-SPIRIT
                    % =====================================================================
                    % initialize arrays
                    mag_NRMSE_SPIRIT_per_iter = zeros(1,S.PARAMS.nIter);  % magnitude recon error
                    dT_RMSE_SPIRIT_per_iter = zeros(1,S.PARAMS.nIter);   % dT (Temperature change) recon error [degree Celsius]
                    dT_NRMSE_SPIRIT_per_iter = zeros(1,S.PARAMS.nIter);   % dT (Temperature change) recon Noramlized error [no units]
                    
                    tic
                    
                    x = x0;
                    x_old = x0;
                    for n=1:S.PARAMS.nIter
                        x = (x + S.GOP*x ).*(mask) + data; % Apply (G-I)*x + x
                        
                        % apply wavelet thresholding
                        X = ifft2c(x); % goto image domain
                        X= zpad(X,ss,ss,nc); % zero-pad to the closest diadic
                        X = W*(X); % apply wavelet
                        X = softThresh(X,wavWeight_SPIRIT); % threshold ( joint sparsity)
                        X = W'*(X); % get back the image
                        X = crop(X,sx,sy,nc); % return to the original size
                        xx = fft2c(X); % go back to k-space
                        x = xx.*mask + data; % fix the data

                        %     if show
                        %         X = ifft2c(x);
                        %         Xsqr = sqrt(sum(abs(X).^2,3));
                        %         figure(20), imshow(Xsqr,[],'InitialMagnification',400);
                        %         title(['TED iter #',num2str(n),' rec']) % 
                        %         drawnow
                        %     end
                        
                       % ============= calc rec errors  ============== 
                        if calc_error_per_iter_flag == 1 | n == S.PARAMS.nIter
                            X_rec = ifft2c(x);
                            Xsqr_rec = sqrt(sum(abs(X_rec).^2,3));
                            [Xsqr_rec] = mag_corrections_per_vendor(Xsqr_rec,S.PARAMS);
                            
                            % ============= magnitude NRMSE  ==============
                            [mag_NRMSE_SPIRIT_per_iter(n),~] = calc_error(Xsqr_gold,Xsqr_rec);

                            % ========== Temperature RMSE ================================
                            % SPIRIT reconstruction
                            [dT_SPIRIT,dT_zoomed_SPIRIT] = TempCalc_v2(DATA_baseline*S.PARAMS.scale_fctr,x*S.PARAMS.scale_fctr,S);
                            dT_zoomed_SPIRIT_err = dT_zoomed_SPIRIT - S.dT_gold_zoomed;
                            [dT_NRMSE_SPIRIT,dT_RMSE_SPIRIT] = calc_error(S.dT_gold_zoomed(:),dT_zoomed_SPIRIT(:));
                            
                            dT_RMSE_SPIRIT_per_iter(n) = dT_RMSE_SPIRIT;
                            dT_NRMSE_SPIRIT_per_iter(n) = dT_NRMSE_SPIRIT;
                            
                            
                            if n==1 | mod(n,10)==0 | n==S.PARAMS.nIter
                                disp(['SPIRIT iter #',num2str(n),' dT RMSE=',num2str(dT_RMSE_SPIRIT)])
                            end
                            
                            
                        elseif n==1 | mod(n,10)==0
                            disp(['SPIRIT iter #',num2str(n)])
                        end
                        
                    end % for n=1:S.PARAMS.nIter
                                    
                    % results to return
                    dT_rec = dT_SPIRIT;
                    dT_rec_zoomed = dT_zoomed_SPIRIT;
                    dT_NRMSE_per_iter = dT_NRMSE_SPIRIT_per_iter;
                    dT_RMSE_per_iter = dT_RMSE_SPIRIT_per_iter;
                    mag_NRMSE_per_iter = mag_NRMSE_SPIRIT_per_iter;
                    
                case 'TED'
                    wavWeight_TED = wavWeight;
                    % =====================================================================
                    %                           TED
                    % =====================================================================
                    
                    % initialize arrays
                    mag_NRMSE_TED_per_iter = zeros(1,S.PARAMS.nIter);     % magnitude recon error
                    dT_RMSE_TED_per_iter = zeros(1,S.PARAMS.nIter);      % dT (Temperature change) recon error [degree Celsius]
                    dT_NRMSE_TED_per_iter = zeros(1,S.PARAMS.nIter);      % dT (Temperature change) recon Noramlized error [no units]
                    
                    x = x0;
                    x_old = x0;
                    X_basline = ifft2c(DATA_baseline); % goto image domain
                    
                    for n=1:S.PARAMS.nIter
                        x = (x + S.GOP*x ).*(mask) + data; % Apply (G-I)*x + x
                        
                        % apply wavelet thresholding
                        X = ifft2c(x); % goto image domain
                        X= zpad(X,ss,ss,nc); % zero-pad to the closest diadic
                        
                        dX = X - X_basline; % calc complex differences
                                                
                        dX_W = W*(dX); % transform the complex differences to the wavelet domain
                                                
                        dX_W_th = softThresh(dX_W,wavWeight_TED); % apply soft thresholding 
                        dX = W'*(dX_W_th); % transform the complex differences back to the image domain
                                                
                        X = X_basline + dX; % compute the complex data for current time frame
                        
                        X = crop(X,sx,sy,nc); % return to the original size
                        xx = fft2c(X); % go back to k-space
                        x = xx.*mask + data; % fix the data
                        
                        %     if show
                        %         X = ifft2c(x);
                        %         Xsqr = sqrt(sum(abs(X).^2,3));
                        %         figure(20), imshow(Xsqr,[],'InitialMagnification',400);
                        %         title(['TED iter #',num2str(n),' rec']) 
                        %         drawnow
                        %     end
                        
                         % ============= calc rec errors  ============== 
                        if calc_error_per_iter_flag == 1 | n == S.PARAMS.nIter
                            
                            X_rec = ifft2c(x);
                            Xsqr_rec = sqrt(sum(abs(X_rec).^2,3));
                            [Xsqr_rec] = mag_corrections_per_vendor(Xsqr_rec,S.PARAMS);
                            
                            % ============= Magnitude RMSE ==============
                            [mag_NRMSE_TED_per_iter(n),~] = calc_error(Xsqr_gold,Xsqr_rec);
                         
                            % ========== Temperature RMSE ================================
                            % TED reconstruction
                            [dT_TED,dT_zoomed_TED] = TempCalc_v2(DATA_baseline*S.PARAMS.scale_fctr,x*S.PARAMS.scale_fctr,S);
                            dT_zoomed_TED_err = dT_zoomed_TED - S.dT_gold_zoomed;
                            [dT_NRMSE_TED,dT_RMSE_TED] = calc_error(S.dT_gold_zoomed(:),dT_zoomed_TED(:));
                            
                            dT_RMSE_TED_per_iter(n) = dT_RMSE_TED;
                            dT_NRMSE_TED_per_iter(n) = dT_NRMSE_TED;
                            
                            if n==1 | mod(n,10)==0 | n==S.PARAMS.nIter
                                disp(['TED iter #',num2str(n),' dT RMSE=',num2str(dT_RMSE_TED)])
                            end
                            
                        elseif n==1 | mod(n,10)==0 | n==S.PARAMS.nIter
                            disp(['TED iter #',num2str(n)])
                        end
                        
                     
                    end % for n=1:S.PARAMS.nIter
                     
                    % results to return
                    dT_rec = dT_TED;
                    dT_rec_zoomed = dT_zoomed_TED;
                    dT_NRMSE_per_iter = dT_NRMSE_TED_per_iter;
                    dT_RMSE_per_iter = dT_RMSE_TED_per_iter; % this is the important measure for Temperature-change images
                    mag_NRMSE_per_iter = mag_NRMSE_TED_per_iter; % this is the imporatant error measure for magnitude images
                    
            end  % switch method_flag
            
            X_rec = ifft2c(x);
            Xsqr_rec = sqrt(sum(abs(X_rec).^2,3));
            [Xsqr_rec] = mag_corrections_per_vendor(Xsqr_rec,S.PARAMS);
            
            
            
            % -------------------------------------------------------
            function x = softThresh(y,t)
                % Soft threshold: S_th(y) = (y/|y|)*max(0,|y|-th)
                % apply joint sparsity soft-thresholding
                % (c) Michael Lustig 2007
                absy = sqrt(sum(abs(y).^2,3));
                
                % compute y/|y|
                unity = y./(repmat(absy,[1,1,size(y,3)])+eps);
                
                % compute max(0,abs(y)-threshold)
                res = absy-t;
                res = (res + abs(res))/2;
                
                % compute (y/|y|)*max(0,|y|-th)
                x = unity.*repmat(res,[1,1,size(y,3)]);
            end % function x = softThresh(y,t)
            
        end % function TED_and_SPIRiT
        
        % ===============================================================================
        
        function S = display_all_slices(S) % data is without scaling!
            % #########################  find useful slices #############################################
            % display gold standard temperature changes - without scaling by scale_fctr
            for slice_i = 1:length(S.PARAMS.slices_for_recon)
                slice = S.PARAMS.slices_for_recon(slice_i);
                FullKspace = squeeze(S.FullKspace_all_slices(:,:,:,:,slice));
                %figure
                MAT = zeros(S.PARAMS.N,S.PARAMS.N,12);  % there's room for 12 images of 256x256
                
                for t_jjj = 1:length(S.PARAMS.t_rec_vec)
                    
                    t_ind = S.PARAMS.t_rec_vec(t_jjj);
                    
                    
                    DATA_baseline = squeeze(FullKspace(:,:,S.PARAMS.t_baseline,:));
                    DATA = squeeze(FullKspace(:,:,t_ind,:));
                    scale_tmp = 1; %10^4;
                    
                    [dT_gold_no_scaling] = TempChangeCalc(DATA_baseline*scale_tmp,DATA*scale_tmp,S.PARAMS);
                    [dT_gold_no_scaling] = dT_corrections_per_vendor(dT_gold_no_scaling,S.PARAMS);
                    
                    % Visual check - for debugging only
                    % figure; imagesc(dT_gold_no_scaling.*S.PARAMS.HIFU_MASK); colormap jet; caxis([-30 30])
                    %   figure; imagesc(dT_gold_no_scaling); colormap jet; caxis([-30 30])
                    
                    MAT(:,:,t_jjj) = dT_gold_no_scaling;
                    %         %figure
                    %         if t_jjj<=9
                    %        % subplot(3,3,t_jjj)
                    %         figure
                    %         imagesc(squeeze(dT_gold_no_scaling))
                    %         axis equal; axis tight; axis off;
                    %         caxis([S.PARAMS.cmin S.PARAMS.cmax]);%caxis([-30 30])
                    %         colormap jet
                    %         c = colorbar;
                    %         c.Label.String = '[^oC]';
                    %         title(['t_i_n_d=',num2str(t_ind)])
                    %         suptitle([S.PARAMS.title,' slice #',num2str(slice),' \DeltaT (full-samp)']) % from t=',num2str(t_baseline),' (baseline) to t=',num2str(t_ind)])
                    %
                    %         end
                end % for t_jjj
                
                % ----- plot if t_jj<=12 ---------
                %                 MAT_row1 = [MAT(:,:,1) MAT(:,:,2) MAT(:,:,3) MAT(:,:,4)];
                %                 MAT_row2 = [MAT(:,:,5) MAT(:,:,6) MAT(:,:,7) MAT(:,:,8)];
                %                 if length(S.PARAMS.t_rec_vec)<=8
                %                     MAT_reshaped = [MAT_row1; MAT_row2;];
                %                 elseif length(S.PARAMS.t_rec_vec)<=12
                %                     MAT_row3 = [MAT(:,:,9) MAT(:,:,10) MAT(:,:,11) MAT(:,:,12)];
                %                     MAT_reshaped = [MAT_row1; MAT_row2; MAT_row3];
                %                 end
% 
%                 % subplot(3,3,t_jjj)
%                 figure
%                 imagesc(squeeze(MAT_reshaped))
%                 axis equal; axis tight; axis off;
%                 caxis([S.PARAMS.cmin S.PARAMS.cmax]);%caxis([-30 30])
%                 colormap jet
%                 c = colorbar;
%                 c.Label.String = '[^oC]';
%                 %title(['t_i_n_d=',num2str(t_ind)])
%                 suptitle([S.PARAMS.title,' slice #',num2str(slice),' \DeltaT (full-samp)']) % from t=',num2str(t_baseline),' (baseline) to t=',num2str(t_ind)])
%                 axis image
                
                cnt = 0;
                for fig_cnt = 1:floor(S.PARAMS.t_rec_vec(end)/12)
                    
                    MAT_tmp = [];
                    for row= 2 %1:3
                        for tt=1:4:12
                            
                            
                            t_jjj = cnt*12+tt;
                            MAT_row = [MAT(:,:,cnt*12+tt) MAT(:,:,cnt*12+tt+1) MAT(:,:,cnt*12+tt+2) MAT(:,:,cnt*12+tt+3)];
                            
                            
                            if (cnt*12+tt+3) == S.PARAMS.t_rec_vec(end)
                                break
                            end
                        end % for tt
                        MAT_tmp = [MAT_tmp; MAT_row];
                    end % for row

%                     figure
%                     imagesc(squeeze(MAT_tmp))
%                     axis equal; axis tight; axis off;
%                     caxis([S.PARAMS.cmin S.PARAMS.cmax]);%caxis([-30 30])
%                     colormap jet
%                     c = colorbar;
%                     c.Label.String = '[^oC]';
%                     %title(['t_i_n_d=',num2str(t_ind)])
%                     suptitle([S.PARAMS.title,' slice #',num2str(slice),' \DeltaT (full-samp)']) % from t=',num2str(t_baseline),' (baseline) to t=',num2str(t_ind)])
%                     axis image
                    
                    cnt = cnt + 1;
                    
                end
                
                
            end
            %end
        end % function S = display_all_slices(S)
       
        
    end % methods
    
end % classdef