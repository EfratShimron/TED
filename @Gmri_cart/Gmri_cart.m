function ob = Gmri_cart(kmask,immask,kfftShift,ImfftShift)

% kmask: mask of sampled k-space locations
if ~islogical(kmask)
    error 'kmask must be a logical nd matrix or 1d vector'
end
ob.kmask = kmask;
if exist('kfftShift','var')
    ob.kfftShift = kfftShift;
else
    ob.kfftShift = true;
end
if exist('ImfftShift','var')
    ob.ImfftShift = ImfftShift;
else
    ob.ImfftShift = true;
end
if exist('immask','var')
    if ~isempty(immask)
        if ~islogical(immask)
            error 'immask must be a logical nd matrix'
        end
        ob.immask = immask;
        if numel(immask) ~= numel(kmask)
            error 'size of immask must = size of kmask'
        end
    end
end
if size(kmask,1) == 1 || size(kmask,2) == 1
  % instead of fftshift, just phase correct after fft (much faster on GPU)
  ob.fftshiftPhs = exp(1i*pi*(0:length(kmask)-1)');
else
  ob.fftshiftPhs = [];
end
ob.adjoint = 0;
ob = class(ob,'Gmri_cart');
