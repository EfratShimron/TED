clear all
close all
clear classes
restoredefaultpath % clear all existing paths - this is important for l1-SPIRIT toolbox!
addpath(genpath(pwd))

% =======================================================================
% -----------------------------------------------------------------------
%                              TED toolbox 
% -----------------------------------------------------------------------
% =======================================================================
% (c) Efrat Shimron, 2020 
% Contact: efrat.s@berkeley.edu
% Follow on Twitter: @Efrat_Shimron

% ------------------------------------------------------------------------
% This is a Matlab toolbox for the Temporal Differences (TED) Compressed
% Sensing methods, published in this paper:
%    E. Shimron, W. Grissom, H. Azhari, "Temporal Differences (TED) Compressed
%    Sensing: A Method for Fast MRgHIFU Temperature Imaging",
%    Accepted for publication NMR in Biomedicine (May 2020)

% This toolbox includes two demos:
% 1. Gel Phantom demo - reproduces Fig. 4 of our paper - data provided by INSIGHTEC, and was acquired with a GE scanner)
% 2. Agar phantom demo - reproduces Fig. 5 of our paper - data is courtsy of Prof. William Grissom,
%                        Vanderbilt University, and was acquired with a Philips Ingenia scanner).

% TED is a general method of dynamic MR-guided-HIFU with parallel imaging.
% Assuming that k-space is fully-sampled at time t=0, TED can accurately reconstruct 
% complex MRI data from sub-sampled k-space data acquired at later time frames (t=1,2...) 

% TED is demonstrated here for temperature reconstruction from subsampled
% data. However, since this method is general, it can be implemented for other types of
% dynamic MRI data. TED is compatible with both single-coil and multi-coil
% data.

% If you publish an interesting implementation - let us know! :-)

% ------------------------------- Acknowledgments -------------------------
% This toolbox was built upon Michael Lustig's l1-SPIRiT toolbox, which was
% downloaded from his website:
% http://people.eecs.berkeley.edu/~mlustig/Software.html
% hence the l1-SPIRiT implementation is included here.
% If you use that method, it is advised to cite one the l1-SPIRiT publication:
% Murphy M, Alley M, Demmel J, Keutzer K, Vasanawala S, Lustig M. "Fast l?-SPIRiT 
% compressed sensing parallel imaging MRI: scalable parallel implementation and 
% clinically feasible runtime". IEEE Transactions on Medical Imaging. 2012;31(6):1250�62.

% This toolbox also includes an implementation of the K-space Hybrid
% Method, which was downloaded from Will Grissom's website:
% https://vuiis.vumc.org/~grissowa/software.html
% If you use it, please cite the following paper:
% Gaur P, Grissom WA. "Accelerated MRI thermometry by direct estimation of 
% temperature from undersampled k-space data". Magnetic Resonance in Medicine. 2015;73(5):1914�1925.



% #################################################################
%                  Comments
% #################################################################
% 1. computing the temperature and the reconstruction error in each
%    iteration really slows down the code, so the code can be made faster
%    by setting calc_error_per_iter_flag=0
% 2. The runtime can be shortened by running the code on GPUs, and using
% "parfor" wherever there is a for loop, especially inside the functions
% fftc and ifftc!

% =============== TO DO - before publicizing the toolbox ==========
% 2. remove the following fields from the function set_params_2020_v2
% (they need to be deleted for every demo separately):
% PARAMS.wavWeight_SPIRIT_vec
% lam_vec
% 3. remove l1-SPIRIT and k-space hybrid functions and their params
% 4. remove the flags SPIRIT_flag, kspace_hybrid_flag, show_flag



demo = 'Gel_phantom_demo'; % Insightec's phantom data
%demo = 'Agar_phantom_demo'; % Grissom's data


SPIRIT_flag = 1;             % set to 1 if you want to run the l1-SPIRiT method
kspace_hybrid_flag = 1;      % set to 1 if you want to run the K-space Hybrid method [Gaur and Grissom, 2015]

show_flag = 0;

calc_error_per_iter_flag = 1;


% set the subsampling reduction factors

switch demo
    case 'Gel_phantom_demo'  % Insightec's phantom
        R_vec = [2:2:10]; % sub-sampling factor
        
    case 'Agar_phantom_demo' % Grissom's data
        R_vec = [2 4 6]; % sub-sampling factor        
  
end


% ==================== Initialize MRgHIFU object ==========================

S = MRgHIFU(demo,R_vec);  % initialize an object of class "MRgHIFU_Series"


S.FLAGS.dT_positivity = 1;  
% This flag is used to enable for a fair comparison with the k-space
% hybrid method, which can reconstruct only non-negative temperature
% changes. In TED, we first reconstruct the temp change (dT), and if this
% flag equals 1 we then null all pixels with a negative phase change. 
% However, this step is not mandatory for the TED method, 



% =========================== load fully sampled data =========================

S = load_kspace_data(S);

% % #########################  display temperature change for all slices #############################################
S = display_all_slices(S); % data is without scaling! 


% %%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
num_slices = length(S.PARAMS.slices_for_recon);

%dT_gold_MAT = zeros(S.PARAMS.N,S.PARAMS.N,length(R_vec),S.PARAMS.NS,S.PARAMS.NT);
%dT_gold_zoomed_MAT = zeros(S.PARAMS.N,S.PARAMS.N,length(R_vec),S.PARAMS.NS,S.PARAMS.NT);

for slice_i =  1:num_slices
    
    slice = S.PARAMS.slices_for_recon(slice_i);
    S.FullKspace = squeeze(S.FullKspace_all_slices(:,:,:,:,slice));
        
        for R_i = 1:length(R_vec)
            R = R_vec(R_i);
            disp(['====================== R=',num2str(R),' ==================='])
            
            S = load_sampling_mask(S,R);
            
            S = scale_and_calc_GOP(S);
            
            % ======================= loop over t_ind ============================
            for t_jjj = 1:length(S.PARAMS.t_rec_vec)
                t_ind = S.PARAMS.t_rec_vec(t_jjj);
                
                disp(['================ t_jjj =',num2str(t_jjj),' ==================='])
                
                % ===================== prepare k-space data & SPIRIT calibration ===================================
                FullKsapce_t_baseline = squeeze(S.FullKspace(:,:,S.PARAMS.t_baseline,:));
                FullKsapce_t_HIFU = squeeze(S.FullKspace(:,:,t_ind,:));
                
                % =================== sample k-space  =================
                SampledKspace = FullKsapce_t_HIFU.*repmat(S.sampling_mask_with_ACS,[1,1,S.PARAMS.NC]); % multiply with sampling matrix
                
                % =================== scale k-space  =================
                % scaling is very useful, since it allows you to
                % reconstruct different datasets with similar compressed
                % sensing parameters. 
                % However, if you use TED (or l1-SPIRiT) for temperature recon,
                % keep in mindt aht since temperature is linearly related to the 
                % phase change, the data must be scaled back to its
                % original scale before the temperature change is computed.
                DATA_baseline = FullKsapce_t_baseline/S.PARAMS.scale_fctr;  % fully sampled k-space data at t=0, scaled
                DATA_full_samp = FullKsapce_t_HIFU/S.PARAMS.scale_fctr; % fully sampled k-space data at t_dynamic, scaled
                DATA_sampled = SampledKspace/S.PARAMS.scale_fctr;  % sub-sampled k-space data at t_dynamic, scaled
                
                % ============== compute dT gold =========================
                [S.dT_gold,S.dT_gold_zoomed] = TempCalc_v2(FullKsapce_t_baseline,FullKsapce_t_HIFU,S);
                dT_gold_MAT(:,:,R_i,slice_i,t_jjj) = S.dT_gold;
                dT_gold_zoomed_MAT(:,:,R_i,slice_i,t_jjj) = S.dT_gold_zoomed;
                
                
                
                
                % =====================================================================
                %                           TED
                % =====================================================================
                method_flag = 'TED';
                wavWeight_TED = S.PARAMS.wavWeight_TED_vec(1);
                [dT_TED,dT_zoomed_TED,dT_NRMSE_TED_per_iter,dT_RMSE_TED_per_iter,mag_NRMSE_TED_per_iter,Xsqr_gold,Xsqr_rec_TED] = TED_and_SPIRiT(DATA_sampled, S, wavWeight_TED,show_flag,DATA_full_samp,DATA_baseline,R,calc_error_per_iter_flag,method_flag);
                
                dT_RMSE_TED_MAT(R_i,t_jjj,slice_i) = dT_RMSE_TED_per_iter(end);
                dT_TED_MAT(:,:,R_i,t_jjj,slice_i) = dT_TED;
                dT_zoomed_TED_MAT(:,:,R_i,t_jjj,slice_i) = dT_zoomed_TED;
                
            
                % =====================================================================
                %                          l1-SPIRIT
                % =====================================================================
                
                if SPIRIT_flag == 1
                    method_flag = 'SPIRIT';
                    wavWeight_SPIRIT = S.PARAMS.wavWeight_SPIRIT_vec(1);
                    [dT_SPIRIT,dT_zoomed_SPIRIT,dT_NRMSE_SPIRIT_per_iter,dT_RMSE_SPIRIT_per_iter,mag_NRMSE_SPIRIT_per_iter,Xsqr_gold,Xsqr_rec_SPIRIT] = TED_and_SPIRiT(DATA_sampled, S, wavWeight_SPIRIT,show_flag,DATA_full_samp,DATA_baseline,R,calc_error_per_iter_flag,method_flag);
                    
                    dT_SPIRIT_MAT(:,:,R_i,t_jjj,slice_i) = dT_SPIRIT;
                    dT_zoomed_SPIRIT_MAT(:,:,R_i,t_jjj,slice_i) = dT_zoomed_SPIRIT;
                    dT_RMSE_SPIRIT_MAT(R_i,t_jjj,slice_i) = dT_RMSE_SPIRIT_per_iter(end);
                    
                end
                  
                
                
                % =====================================================================
                %                           K-space Hybrid Method
                % =====================================================================
                if kspace_hybrid_flag == 1
                    
               
                        [dT_zoomed_Gaur] = Hybrid(S.FullKspace,S.PARAMS,S.sampling_mask_with_ACS,t_jjj,S.PARAMS.lam);
                        
                        
                        dT_zoomed_Gaur_err = dT_zoomed_Gaur - S.dT_gold_zoomed;
                        [NRMSE_Gaur,RMSE_Gaur] = calc_error(S.dT_gold_zoomed(:),dT_zoomed_Gaur(:));
                        
                        dT_RMSE_Gaur_MAT(R_i,t_jjj,slice_i) = RMSE_Gaur;
                        dT_zoomed_Gaur_MAT(:,:,R_i,t_jjj,slice_i) = dT_zoomed_Gaur;
                    
                else
                    dT_zoomed_Gaur = -20*ones(size(S.dT_gold_zoomed));
                    dT_zoomed_Gaur_err = -20*ones(size(S.dT_gold_zoomed));
                    RMSE_Gaur = 0;
                    
                end % if kspace_hybrid_flag
                
                
                
                 
                             
                
            end % for t_jjj
        end % for R_i
        
    
end % for slice_i


% % ##################################################################
% %                    plots 
% % ##################################################################

switch demo
    case 'Gel_phantom_demo'
        text_shift = 118;  % horixontal shift for showing the RMSE values on the figure
        
    case 'Agar_phantom_demo'
       text_shift = 67; %V*4+10;
end

if SPIRIT_flag==1 & kspace_hybrid_flag==1
    
    SPACE = -20*ones(length(S.PARAMS.x1:S.PARAMS.x2),1);
    
    for slice_i =  1:num_slices
        slice = S.PARAMS.slices_for_recon(slice_i);
        
        for t_jjj = 1:length(S.PARAMS.t_rec_vec)
            t_ind = S.PARAMS.t_rec_vec(t_jjj);
            MAT = [];
            ERR_MAT = [];
            
            dT_gold_zoomed = dT_gold_zoomed_MAT(:,:,R_i,t_jjj,slice_i);   % squeeze(dT_gold_zoomed_MAT(:,:,t_jjj));  % xxx FIX THIS
            for R_i = 1:length(R_vec)
                dT_zoomed_Gaur = squeeze(dT_zoomed_Gaur_MAT(:,:,R_i,t_jjj,slice_i));
                dT_zoomed_SPIRIT = squeeze(dT_zoomed_SPIRIT_MAT(:,:,R_i,t_jjj,slice_i));
                dT_zoomed_TED = squeeze(dT_zoomed_TED_MAT(:,:,R_i,t_jjj,slice_i));
                
                dT_zoomed_Gaur_err = dT_zoomed_Gaur - dT_gold_zoomed;
                dT_zoomed_SPIRIT_err = dT_zoomed_SPIRIT - dT_gold_zoomed;
                dT_zoomed_TED_err = dT_zoomed_TED - dT_gold_zoomed;
                
                if R_i == 1
                    MAT = [dT_gold_zoomed SPACE dT_zoomed_Gaur SPACE dT_zoomed_SPIRIT SPACE dT_zoomed_TED];
                else
                    MAT = [MAT ; -20*ones(size(dT_gold_zoomed)) SPACE dT_zoomed_Gaur SPACE dT_zoomed_SPIRIT SPACE dT_zoomed_TED];
                end
                ERR_MAT = [ERR_MAT; dT_zoomed_Gaur_err SPACE dT_zoomed_SPIRIT_err SPACE dT_zoomed_TED_err ];
                
            end
            
            H = size(dT_gold_zoomed,2);  % horizontal shift for each method
            V = size(dT_gold_zoomed,1);  % vertical shift for each method
            
            VERTICAL_SPACE = -20*ones(size(MAT,1),4);
            figure('Name',['RECS:', S.PARAMS.title,' slice #',num2str(slice),' t=',num2str(t_ind)])
            imagesc([MAT VERTICAL_SPACE ERR_MAT])
            axis equal; axis tight; axis off;
            caxis([S.PARAMS.cmin S.PARAMS.cmax]);
            colormap jet; colorbar;
            title('Rerence + recs; errors')
            hold on
            text(5,3,[sprintf('Ref')],'Color','k','FontSize',10)
            hold on
            text(2+H,3,[sprintf('Hybrid')],'Color','k','FontSize',10)
            hold on
            text(2+2*H+1,3,[sprintf('SPIRiT')],'Color','k','FontSize',10)
            hold on
            text(2+3*H+5,3,[sprintf('TED')],'Color','k','FontSize',10)
            hold on
            text(text_shift,3,[sprintf('Hybrid')],'Color','k','FontSize',10)
            hold on
            text(text_shift+H+1,3,[sprintf('SPIRiT')],'Color','k','FontSize',10)
            hold on
            text(text_shift+2*H+5,3,[sprintf('TED')],'Color','k','FontSize',10)
            for R_i = 1:length(R_vec)
                RMSE_Gaur = dT_RMSE_Gaur_MAT(R_i,t_jjj,slice_i);
                RMSE_SPIRIT = dT_RMSE_SPIRIT_MAT(R_i,t_jjj,slice_i);
                RMSE_TED = dT_RMSE_TED_MAT(R_i,t_jjj,slice_i);
                hold on
                text(text_shift+1,R_i*(S.PARAMS.x2-S.PARAMS.x1)-0.5,[sprintf('%.3f',RMSE_Gaur)],'Color','k','FontSize',10)
                hold on
                text(text_shift+H,R_i*(S.PARAMS.x2-S.PARAMS.x1)-0.5,[sprintf('%.3f',RMSE_SPIRIT )],'Color','k','FontSize',10)
                hold on
                text(text_shift+2*H,R_i*(S.PARAMS.x2-S.PARAMS.x1)-0.5,[sprintf('%.3f',RMSE_TED )],'Color','k','FontSize',10)
            end
            c = colorbar;
            c.Label.String = '[^oC]';
            suptitle([S.PARAMS.title,' Time Frame #',num2str(t_ind),' slice #',num2str(slice)])
            
        end % for t_jjj
        
    end % for slice_i
end


% % ================== full-FOV phase map =======================
switch demo
    case 'Gel_phantom_demo'
        R_i_wanted = length(R_vec);
        slice4display = 2;
        t_jjj_tmp = find(S.PARAMS.t_example==S.PARAMS.t_rec_vec);
        gold_full_FOV = squeeze(dT_gold_MAT(:,:,R_i_wanted,slice4display,t_jjj_tmp));
        TED_full_FOV = squeeze(dT_TED_MAT(:,:,R_i_wanted,slice4display,t_jjj_tmp));
        
        % gold standard
        figure; 
        imagesc(squeeze(gold_full_FOV))
        axis equal; axis tight; axis off;
        caxis([S.PARAMS.cmin S.PARAMS.cmax]);%caxis([-30 30])
        colormap jet
        c = colorbar;
        c.Label.String = '[^oC]';
        title('Reference (gold standard) phase map')
        %title(['Reference', S.PARAMS.title,', Full-FOV Gold Standard dT for t=',num2str(S.PARAMS.t_example),' R=',num2str(R_vec(R_i_wanted)),' slice #',num2str(slice4display)]
        
        % TED
        figure;
        imagesc(squeeze(TED_full_FOV))
        axis equal; axis tight; axis off;
        caxis([S.PARAMS.cmin S.PARAMS.cmax]);%caxis([-30 30])
        colormap jet
        c = colorbar;
        c.Label.String = '[^oC]';
        title('TED phase map')
        %title(['TED ',S.PARAMS.title,', Full-FOV TED dT for t=',num2str(S.PARAMS.t_example),' R=',num2str(R_vec(R_i_wanted)),' slice #',num2str(slice4display)]

end


