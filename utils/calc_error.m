function [NRMSE,RMSE] = calc_error(correct_im,rec)
% Compute the Root Mean Squared Error (RMSE) and the Normalized Root Mean Squared Error (NRMSE)
% (c) Efrat Shimron, 2019

% s_0 = abs(correct_im);
% s_rec = abs(rec);

s_0 = correct_im;
s_rec = rec;

% switch flag 
%     case 'positive_vals_only'
%         s_0(s_0<=0) = 0;
%         s_rec(s_rec<0)=0;
% end

% % throw away pixels for which the true image == 0
% nnz_inds = find(s_0~=0);
% s_0 = s_0(nnz_inds);
% s_rec = s_rec(nnz_inds);

len = length(s_0);
RMSE = sqrt(sum((s_rec(:)-s_0(:)).^2 ) / (size(s_0,1)*size(s_0,2)));
dd = ( max(s_0(:)) - min(s_0(:))  );
NRMSE = RMSE/dd;

end