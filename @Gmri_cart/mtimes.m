function res = mtimes(a,bb)

if a.adjoint
    % k -> image, via ifft
    if min(size(a.kmask)) > 1 
        % all-dim ifft
        kmat = zeros(size(a.kmask),class(bb));
        kmat(a.kmask) = bb;
        if a.kfftShift == true
            kmat = fftshift(kmat);
        end
        res = ifftn(kmat);
        if a.ImfftShift == true
          res = fftshift(res);
        end
    end
    if size(a.kmask,1) == 1
        % only apply ifft in column dimension
        kmat = zeros(length(a.kmask),class(bb));
        kmat(:,a.kmask) = reshape(bb,[length(a.kmask) sum(a.kmask)]);
        if a.kfftShift == true
            kmat = fftshift(kmat,2);
        end
        if a.ImfftShift == true
          kmat = bsxfun(@times,kmat,conj(a.fftshiftPhs).');
        end
        res = ifft(kmat,[],2);
        %if a.ImfftShift == true
        %  res = fftshift(res,2);
        %end
    end
    if size(a.kmask,2) == 1
        % only apply ifft in row dimension
        kmat = zeros(length(a.kmask),class(bb));
        kmat(a.kmask,:) = reshape(bb,[sum(a.kmask) length(a.kmask)]);
        if a.kfftShift == true
            kmat = fftshift(kmat,1);
        end
        if a.ImfftShift == true
          kmat = bsxfun(@times,kmat,conj(a.fftshiftPhs));
        end
        res = ifft(kmat,[],1);
        %if a.ImfftShift == true
        %  res = fftshift(res,1);
        %end
    end
    if isfield(a,'immask');
        res = res(a.immask);
    else
        res = res(:);
    end
else
    % image -> k, via fft
    if isfield(a,'immask')
        tmp = zeros(numel(a.immask),1,class(bb));
        tmp(a.immask) = bb;
        bb = tmp;
    end
    if min(size(a.kmask)) > 1 
        % all-dim fft
        if a.ImfftShift == true
          res = fftn(fftshift(reshape(bb,size(a.kmask))))./numel(a.kmask);
        else
          res = fftn(reshape(bb,size(a.kmask)))./numel(a.kmask);
        end
        if a.kfftShift == true
            res = fftshift(res);
        end
        res = res(a.kmask);
    end
    if size(a.kmask,1) == 1
        % only apply fft in column dimension
        if a.ImfftShift == true
          res = bsxfun(@times,fft(reshape(bb,[length(a.kmask) length(a.kmask)]),[],2),a.fftshiftPhs.')./length(a.kmask);
          %res = fft(fftshift(reshape(bb,[length(a.kmask) length(a.kmask)]),2),[],2)./length(a.kmask);
        else
          res = fft(reshape(bb,[length(a.kmask) length(a.kmask)]),[],2)./length(a.kmask);
        end
        if a.kfftShift == true
            res = fftshift(res,2);
        end
        res = res(:,a.kmask);
        res = res(:);
    end
    if size(a.kmask,2) == 1
        % only apply fft in row dimension
        if a.ImfftShift == true
          res = bsxfun(@times,fft(reshape(bb,[length(a.kmask) length(a.kmask)]),[],1),a.fftshiftPhs)./length(a.kmask);
          %res = fft(fftshift(reshape(bb,[length(a.kmask) length(a.kmask)]),1),[],1)./length(a.kmask);
        else
          res = fft(reshape(bb,[length(a.kmask) length(a.kmask)]),[],1)./length(a.kmask);
        end
        if a.kfftShift == true
            res = fftshift(res,1);
        end
        res = res(a.kmask,:);
        res = res(:);
    end
end

