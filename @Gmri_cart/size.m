function sz = size(a)

if min(size(a.kmask)) > 1
    sz = [sum(a.kmask(:)) numel(a.kmask)];
else
    sz = [sum(a.kmask)*length(a.kmask) length(a.kmask)*length(a.kmask)]; 
end
