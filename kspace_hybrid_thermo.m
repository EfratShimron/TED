function [theta,A,c,f,Ac,algp] = kspace_hybrid_thermo(acqp,thetainit,algp)

%|function kspace_hybrid_thermo
%|
%| Inputs:
%|  acqp    Acquisition parameters structure containing (required):
%|              data        [Nk,Nc]       Nc complex k-space data vectors
%|              k           [Nk,Nd]       Nd k-space sample vectors (cycles/cm)
%|                       OR [Nkx,Nky,Nkz] logical Cartesian k-space sampling mask
%|              fov         1             Field of view (cm) (Non-Cartesian only)
%|              mask        [Nx,Ny,Nz]    Binary mask over the FOV (Non-Cartesian only)
%|              L           [Nx*Ny*Nz*Nc,Nl] Multibaseline image library
%|  thetainit   [Nx,Ny,Nz]  Initial temperature map (real, negative). Real part is
%|                          temperature-induced phase at TE. (optional)
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
%| Outputs:
%|  theta       [Nx,Ny,Nz]    Complex temperature map
%|  A           [Nx*Ny*Nz,Np] Polynomial matrix (may be masked)
%|  c           [Np,Nc]       Polynomial coeffs
%|  f           [Nx,Ny,Nz,Nc] Baseline estimate
%|  Ac          [Nx,Ny,Nc]    Polynomial phase estimate (embedded into original mask)
%|  algp        struct        Final algorithm parameters structure
%|
%| Copyright 2014-05-19, William A Grissom, Pooja Gaur, Vanderbilt University

% ----------------- version 2: --------------------
% theta is not required to be only positive, i.e.
% the line 
% theta(theta >= 0) = 0; is disabled


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('algp','var')
    algp = struct();
end
if ~isfield(algp,'order')
    algp.order = 0; % zeroth-order only (phase drift)
end
if ~isfield(algp,'lam')
    algp.lam = [10^-6 10^-6]; % very small values
    if algp.modeltest
    	algp.lam = [0 10^6];
    end
end
if ~isfield(algp,'beta')
    algp.beta = -1; % turn off roughness penalty if beta not supplied
end
if ~isfield(algp,'dofigs')
    algp.dofigs = 0;
end
if ~isfield(algp,'thiters')
    algp.thiters = 10; % theta iterations
end
if ~isfield(algp,'citers')
    algp.citers = 5; % c iterations
end
if ~isfield(algp,'bls')
    algp.bls = 1; % use boyd's backtracking line search
end
if ~isfield(algp,'masknz')
    algp.masknz = [];
else
    if ~isempty(algp.masknz)
        algp.masknz = algp.masknz(acqp.mask);
    end
end
if ~isfield(algp,'modeltest')
    algp.modeltest = 0;
end
if ~isfield(algp,'sumMask')
    algp.sumMask = false;
end
if ~isfield(algp,'stopFrac')
    algp.stopFrac = 0.001;
end
if ~isfield(algp,'useGPU')
  algp.useGPU = false;
end
if ~isfield(acqp,'dcf')
    acqp.dcf = 1; % optional density compensation function
end

disp('Performing k-space hybrid thermometry.');

Nc = size(acqp.data,2); % Number of rx coils

if islogical(acqp.k)
    disp('k-space is logical array; using Gmri_cart.');
    acqp.fov = 0;
    % pre-fftshift the k-data if Cartesian, so we don't have to do it
    % iteratively
    if min(size(acqp.k)) > 1
        acqp.mask = true(size(acqp.k));
        if ~exist('thetainit','var')
            thetainit = zeros(size(acqp.k));
        end
        kShift = fftshift(acqp.k);
        for ii = 1:Nc
            dataShift = zeros(size(acqp.k));
            dataShift(acqp.k) = acqp.data(:,ii);
            dataShift = fftshift(dataShift);
            acqp.data(:,ii) = dataShift(kShift);
            clear dataShift
        end
        acqp.k = kShift; % Now that the data is shifted, we need to always use
        clear kShift
        % shifted k mask
    end
    if size(acqp.k,1) == 1 % we have undersampling only in column dim
        acqp.mask = true(length(acqp.k)); % assume square k-space and images for Cartesian
        if ~exist('thetainit','var')
            thetainit = zeros(length(acqp.k));
        end
        % fft the data to image domain in row dim and
        % pre-fftshift the k-space data in column dim
        kShift = fftshift(acqp.k);
        for ii = 1:Nc
            dataShiftHybrid = zeros(length(acqp.k)); % assume square k-space + images
            dataShiftHybrid(:,acqp.k) = reshape(acqp.data(:,ii),[length(acqp.k) sum(acqp.k)]);
            dataShiftHybrid = fftshift(ifft(fftshift(dataShiftHybrid),[],1),1)*length(acqp.k); % ifft in row dim
            acqp.data(:,ii) = col(dataShiftHybrid(:,kShift));
            clear dataShiftHybrid;
        end
        % fftshift the mask
        acqp.k = kShift;
        clear kShift
    end
    if size(acqp.k,2) == 1 % we have undersampling only in row dim
        acqp.mask = true(length(acqp.k)); % assume square k-space and images for Cartesian
        if ~exist('thetainit','var')
            thetainit = zeros(length(acqp.k));
        end
        % fft the data to image domain in column dim and
        % pre-fftshift the k-space data in row dim
        kShift = fftshift(acqp.k);
        for ii = 1:Nc
            dataShiftHybrid = zeros(length(acqp.k)); % assume square k-space + images
            dataShiftHybrid(acqp.k,:) = reshape(acqp.data(:,ii),[sum(acqp.k) length(acqp.k)]);
            dataShiftHybrid = fftshift(ifft(fftshift(dataShiftHybrid,2),[],2),2)*length(acqp.k); % ifft in column dim
            dataShiftHybrid = fftshift(dataShiftHybrid,1);
            acqp.data(:,ii) = col(dataShiftHybrid(kShift,:));
            clear dataShiftHybrid;
        end
        % fftshift the mask
        acqp.k = kShift;
        clear kShift
    end
else
    disp('k-space is double array; using Gmri');
    if ~exist('thetainit','var')
        thetainit = zeros(sum(acqp.mask(:)),1);
    end
end

Ns = sum(acqp.mask(:)); % Number of spatial locations

%%%%%%%%%%%%%%%%%%%%%%%%%
% Build objects
%%%%%%%%%%%%%%%%%%%%%%%%%

% build polynomial matrix
A = buildA(acqp.mask,algp.order);

% build system matrix
G = buildG(Nc,acqp.k,acqp.fov,acqp.mask);

% build penalty object
if algp.beta > 0
    R = Robject(acqp.mask,'order',2,'beta',algp.beta,'type_denom','matlab');
else
    R = [];
end

if algp.useGPU && islogical(acqp.k)
  % use GPU for processing (experimental). Currently only works with 2DFT (i.e. Gmri_cart),
  % and the fftshifts are much slower on GPU.
  acqp.data = gpuArray(acqp.data);
  acqp.L = gpuArray(acqp.L);
  acqp.dcf = gpuArray(acqp.dcf);
  acqp.mask = gpuArray(acqp.mask);
  A = gpuArray(A);
  algp.lam = gpuArray(algp.lam);
  R = gpuArray(R);
  thetainit = gpuArray(thetainit);
elseif algp.useGPU && ~islogical(acqp.k)
  warning 'Ignoring useGPU flag since data is not Cartesian'
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% get initial f,c,theta
%%%%%%%%%%%%%%%%%%%%%%%%%

% get initial baseline estimate f
f = f_update(acqp.data,zeros(Ns,1),zeros(Ns,1),acqp.L,G,acqp.dcf);
c = zeros(size(A,2),1,class(A));
Ac = A*c;

% get initial theta
theta = thetainit(acqp.mask(:));

if algp.modeltest
  theta = zeros(Ns,1,class(theta));
end

% force negativity
% theta(theta >= 0) = 0; % disabled in version 2

% Normalize data to get on scale with penalties
dnorm = median(sqrt(acqp.dcf(:)).*abs(acqp.data(:)))*sqrt(length(acqp.data(:)));
if dnorm ~= 0
  acqp.data = acqp.data / dnorm;
  f = f / dnorm;
  acqp.L = acqp.L / dnorm;
else
  disp ['Warning: normalization = 0, so not applied. This can ' ...
        'happen when the object has been masked. lam ' ...
        'may need tweaking.'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1-Penalized Component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(algp.masknz) % only run if we don't already have a heating mask

    costOld = Inf;
    cost = cost_eval(acqp.data,G,f,Ac,theta,R,algp.lam,acqp.dcf);
    itr = 0;
    fprintf('L1-penalized iteration %d, cost = %f\n',itr,cost);
    while costOld-cost >= algp.stopFrac*costOld

        % update baseline image estimates
        f = f_update(acqp.data,Ac,theta,acqp.L,G,acqp.dcf);

        % update poly coeffs
        c = c_update(acqp.data,A,c,theta,f,G,algp,acqp.dcf);
        Ac = A*c;

        % update temp shift
        if ~algp.modeltest
            theta = theta_update(acqp .data,Ac,theta,f,G,1,algp,algp.lam,R,[],acqp.dcf,algp.modeltest);

            if algp.dofigs;figure(201);
                tmp = zeros(size(acqp.mask));tmp(acqp.mask) = -real(theta);
                subplot(121); imagesc(tmp); axis image; title 'Estimated phase';
                tmp = zeros(size(acqp.mask));tmp(acqp.mask) = -real(theta) > 0;
                subplot(122); imagesc(tmp); axis image; title 'Significant phase'
                drawnow;
            end
        end

        % calculate cost with updated parameters
        costOld = cost;
        cost = cost_eval(acqp.data,G,f,Ac,theta,R,algp.lam,acqp.dcf);

        itr = itr + 1;
        fprintf('L1-penalized iteration %d, cost = %f\n',itr,cost);

    end

    % get a mask of potential temperature shifts.
    algp.masknz = theta ~= 0;

    if algp.modeltest
        tmp = sqrt(sum(abs(reshape(f,[Ns Nc])).^2,2));
        algp.masknz = tmp > 0.1*max(tmp);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masked Component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(algp.masknz)

    % run theta_update for nonzero pixels, with no sparsity regularization.
    % we do this because we know that sparsity regularization will
    % attenuate the map somewhat, so we need to relax that effect.
    costOld = Inf;
    cost = cost_eval(acqp.data,G,f,Ac,theta,R,[0 0],acqp.dcf);
    itr = 0;
    fprintf('Masked iteration %d, cost = %f\n',itr,cost);
    while costOld-cost >= algp.stopFrac*costOld

        if ~algp.modeltest
            % update image estimate
            f = f_update(acqp.data,Ac,theta,acqp.L,G,acqp.dcf);

            % update poly coeffs
            c = c_update(acqp.data,A,c,theta,f,G,algp,acqp.dcf);
            Ac = A*c;
        end

        % update temp shift
        theta = theta_update(acqp.data,Ac,theta,f,G,algp.masknz,algp,[0 0],R,algp.sumMask,acqp.dcf,algp.modeltest);

        if algp.dofigs;figure(201);
            tmp = zeros(size(acqp.mask));tmp(acqp.mask) = -real(theta);
            subplot(121); imagesc(tmp); axis image; title 'Estimated phase';
            tmp = zeros(size(acqp.mask));tmp(acqp.mask) = -real(theta) > 0;
            subplot(122); imagesc(tmp); axis image; title 'Significant phase'
            drawnow;
        end

        % calculate cost with updated parameters
        costOld = cost;
        cost = cost_eval(acqp.data,G,f,Ac,theta,R,[0 0],acqp.dcf);

        itr = itr + 1;
        fprintf('Masked Iteration %d, cost = %f\n',itr,cost);

    end
end

% embed final results into full image matrix
tmp = zeros(size(acqp.mask)); 
tmp(acqp.mask) = gather(theta);
theta = tmp; % theta is the output of this function 
tmp = zeros([size(acqp.mask) Nc]);
tmp(repmat(acqp.mask,[ones(1,length(size(acqp.mask))) Nc])) = gather(f);
f = tmp;
tmp = zeros(size(acqp.mask)); 
tmp(acqp.mask) = gather(Ac);
Ac = tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Build the polynomial matrix
%
function A = buildA(mask,order)

if length(size(mask)) == 2 % build a 2D polynomial matrix

    [yc,xc] = meshgrid(linspace(-1/2,1/2,size(mask,2)), ...
        linspace(-1/2,1/2,size(mask,1)));
    yc = yc(:);
    xc = xc(:);
    A = [];
    for yp = 0:order
        for xp = 0:(order-yp)
            A = [A (xc.^xp).*(yc.^yp)];
        end
    end
    A = A(mask(:),:);

else % build a 3D polynomial matrix

    [zc,yc,xc] = meshgrid(linspace(-1/2,1/2,size(mask,3)), ...
        linspace(-1/2,1/2,size(mask,2)), linspace(-1/2,1/2,size(mask,1)));
    zc = zc(:);
    yc = yc(:);
    xc = xc(:);
    A = [];
    for yp = 0:order
        for xp = 0:(order-yp)
            for zp = 0:(order-(yp+xp))
                A = [A (xc.^xp).*(yc.^yp).^(zc.^zp)];
            end
        end
    end
    A = A(mask(:),:);

end

%
% Build the system matrices
%
function G = buildG(Nc,k,fov,mask)

if ~isempty(k) % image-domain
    if ~islogical(k) % non-cartesian

        % build system matrix
        if size(k,3) == 1 % 1 shot
            G = Gmri(k,mask,'fov',fov,'basis',{'dirac'});
        else % multishot
            nshot = size(k,3);
            for ii = 1:nshot % build a system matrix for each shot
                Gsub{ii} = Gmri(k(:,:,ii),mask,'fov',fov,'basis',{'dirac'});
            end
            G = block_fatrix(Gsub,'type','col');
        end

    else % cartesian

        kfftShift = false; % switch to do second fftshift, in frequency domain
        % If false, must make sure k-space data is not centered before starting
        % algorithm
        G = Gmri_cart(k,[],kfftShift);

    end

    if Nc > 1 % multiple coils; replicate the nufft's into a block-diag matrix
        for ii = 1:Nc
            tmp{ii} = G;
        end
        G = block_fatrix(tmp,'type','diag');
    end
else
    G = 1; % image domain - no transform
end


%
% Evaluate cost
%
function cost = cost_eval(data,G,f,Ac,theta,R,lam,dcf)

err = data(:) - G*col(bsxfun(@times,f,exp(1i*(Ac+theta))));
cost = 1/2*real(err'*(dcf(:).*err));

if exist('R','var')
    if ~isempty(R)
        cost = cost + R.penal(R,real(theta)) + R.penal(R,imag(theta));
    end
end
if exist('lam','var')
    if ~isempty(lam)
        cost = cost - lam(1)*sum(real(theta));
        if lam(2) > 0
            cost = cost + lam(2)*sum(imag(theta));
        end
    end
end

%
% Update heat phase shift vector theta
%
function theta = theta_update(data,Ac,theta,f,G,masknz,algp,lam,R,sumMask,dcf,modeltest)

% Polak-Ribiere PCG algorithm from JA Fessler's book, chapter 2, 11.7.13
g = [];
thresh = pi/1000;
for nn = 1:algp.thiters
    gold = g;
    g = masknz.*gradcalc_theta(data,Ac,theta,f,G,R,lam,dcf);
    if exist('sumMask','var')
        % update phase equally over hot spot during masked iterations,
        % since we expect the l1 penalty to just shrink everything by the
        % same amount
        if sumMask == true
            g = sum(g(:))*masknz;
        end
    end
    if nn == 1
        dir = -g;
    else
        gamma = max(0,real(g'*(g-gold))/real(gold'*gold));
        dir = -g + gamma*dir;
    end
    % dir(dir > 0 & -theta < thresh) = 0; % Fessler 11.11.1
    % WAG 10-18-2016: when running modeltest, also set lam(1) = 0; lam(2) = 10^6
    if ~modeltest
      dir(real(dir) > 0 & -real(theta(:,end)) < thresh) = 1i*imag(dir(real(dir) > 0 & -real(theta(:,end)) < thresh));
      dir(imag(dir) < 0 & imag(theta(:,end)) < thresh) = real(dir(imag(dir) < 0 & imag(theta(:,end)) < thresh));
    end

    [t,breakOut] = stepcalc_theta(dir,data,Ac,theta,f,G,R,lam,100*min(1,pi/2/max(abs(dir))),algp,g,dcf);

    if ~modeltest
      z = theta + t*dir;
      if any(real(z) > 0) || any(imag(z) < 0)
          %dir = z.*(z < 0) - theta;
          dir = real(z).*(real(z) < 0) - real(theta(:,end)) + 1i*(imag(z).*(imag(z) > 0) - imag(theta(:,end)));
          [t,breakOut] = stepcalc_theta(dir,data,Ac,theta,f,G,R,lam,1,algp,g,dcf);
      end
    end
    if breakOut == true;break;end
    theta = theta + t*dir;
end

%
% Calculate gradient of cost wrt theta
%
function g = gradcalc_theta(data,Ac,theta,f,G,R,lam,dcf)

% data fidelity derivatives
Nc = size(f,2);
img = col(bsxfun(@times,f,exp(1i*(Ac+theta))));
g = 1i*sum(reshape(conj(img).*(G'*(dcf(:).*(data(:) - G*img))),[length(theta) Nc]),2);
if lam(2) <= 0
    g = real(g) - lam(1); % l1 penalty derivatives; real theta
else
    g = real(g) - lam(1) + 1i*(imag(g) + lam(2)); % l1 penalty derivatives;
    % complex theta
end

if ~isempty(R) % roughness penalty derivatives
    g = g + R.cgrad(R,real(theta)) + 1i*R.cgrad(R,imag(theta));
end

%
% Calculate theta step size
%
function [t,breakOut] = stepcalc_theta(dir,data,Ac,theta,f,G,R,lam,tmax,algp,thetagrad,dcf)

% use boyd's backtracking line search, which usually requires fewer cost evaluations

% calculate current cost
cost = cost_eval(data,G,f,Ac,theta,R,lam,dcf);

% line search to get step
costt = cost;
a = 0.5; b = 0.5; t = tmax/b;
while (costt > cost + a*t*real(thetagrad'*dir)) && t > 10^-6

    % reduce t
    t = b*t;

    % get test point
    thetat = theta + t*dir;

    % calculate cost of test point
    costt = cost_eval(data,G,f,Ac,thetat,R,lam,dcf);

end

if t == tmax/b % loop was never entered; return zero step
    t = 0;
end
if cost-costt >= algp.stopFrac*cost
    breakOut = false;
else
    breakOut = true;
end


%
% Update polynomial coefficient vector c
%
function c = c_update(data,A,c,theta,f,G,algp,dcf)

g = []; % gradient
for nn = 1:algp.citers
    gold = g;
    g = gradcalc_c(data,A,c,theta,f,G,dcf);
    if nn == 1
        dir = -g;
    else
        gamma = max(0,real(g'*(g-gold))/real(gold'*gold));
        dir = -g + gamma*dir;
    end
    alpha = stepcalc_c(dir,data,A,c,theta,f,G,min(1,pi/2/max(abs(dir(:)))),algp.bls,g,dcf);
    c = c + alpha*dir;
end


%
% Calculate gradient of cost wrt c
%
function g = gradcalc_c(data,A,c,theta,f,G,dcf)

Nc = size(f,2);
img = col(bsxfun(@times,f,exp(1i*(A*c+theta))));
g = A'*real(sum(reshape(1i*conj(img).*(G'*(dcf(:).*(data(:) - G*img))),...
    [length(theta) Nc]),2));


%
% Calculate step size for c
%
function t = stepcalc_c(dir,data,A,c,theta,f,G,tmax,bls,cgrad,dcf)

% use boyd's backtracking line search, which usually requires fewer cost evaluations

% calculate current cost
cost = cost_eval(data,G,f,A*c,theta,[],[],dcf);

% line search to get step
costt = cost;
a = 0.5; b = 0.5; t = tmax/b;
while (costt > cost + a*t*real(cgrad'*dir)) && t > 10^-6

  % reduce t
  t = b*t;

  % get test point
  ct = c + t*dir;

  % calculate cost of test point
  costt = cost_eval(data,G,f,A*ct,theta,[],[],dcf);

end

if t == tmax/b % loop was never entered; return zero step
  t = 0;
end

%
% Update baseline image estimates
%
function f = f_update(data,Ac,theta,L,G,dcf)

Nc = size(data,2); % # coils

if size(L,2) > 1 % if more than one library image

    % project library images to k-space
    for ii = 1:size(L,2)
        Lk(:,ii) = sqrt(dcf(:)).*(G*(L(:,ii).*repmat(exp(1i*(Ac+theta)),[Nc 1])));
    end
    LtL = real(Lk'*Lk);

    % set up constraints
    Ceq = ones(1,size(L,2));
    beq = 1;

    % set up cost
    c = -real((data(:).*sqrt(dcf(:)))'*Lk);

    % solve
    options = optimset;options.MaxIter = 100000;options.Algorithm = 'active-set';
    wts = quadprog(double(LtL),double(c),[],[],Ceq,beq,zeros(size(L,2),1),ones(size(L,2),1),[],options);

    % get f
    f = L * wts;

else

    % only one baseline, so weight vector = 1;
    f = reshape(L,[length(L)/Nc Nc]);

end

%
% Helper function
%
function out = col(in)

out = in(:);
