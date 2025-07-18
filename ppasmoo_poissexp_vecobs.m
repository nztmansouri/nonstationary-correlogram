% Point-process adaptive smoothing w/ Poisson likelihood (log-link)
%  filtering via Eden et al. Neural Comp 2004
%  then a backward pass based on Rauch-Tung-Striebel
%  modified to use vector observations

function [predNLL,gamma,Sigma,lam] = ppasmoo_poissexp_vecobs(n,X,b0,W0,F,Q,offset,doSmooth)

% For vector observations of dim k, modeled on T timepoints with p predictors
%
% Inputs...
%   n is [T x k] matrix
%   X is [T x k x p] matrix of corresponding predictors
%   b0 is [p x 1] vector of initial parameters for t=0
%   W0 is [p x p] matrix of initial covariance for t=0
%   F is [p x p] linear process model with noise covariance Q [p x p]
%   (optional) offset is [T x k] matrix for the final poisson model
%   (optional) doSmooth boolean (default true)
%
% Outputs
%   gamma [T x p] mean parameters
%   Sigma [T x p x p] parameter covariance
%   lam [T x k] predictions
%   predLL - predictive log likelihood (for optimizing Q)

if nargin<7 || isempty(offset), offset=zeros(size(n)); end
if nargin<8, doSmooth=true; end

% Preallocate
gamma   = zeros(length(b0),size(n,1));
Sigma   = zeros([size(W0) size(n,1)]);
lam = n*0;

% Initialize
gamma(:,1)   = b0;
Sigma(:,:,1) = W0;
lam(1,:)   = exp((squeeze(X(1,:,:))*b0) + offset(1,:)')';

gamma_pred = gamma;
Sigma_pred = Sigma;

warning('')
fprintf('Filtering...')
tic
% Forward-Pass (Filtering)
for i=2:size(n,1)
    gamma_pred(:,i) = F*gamma(:,i-1);
    Xi = squeeze(X(i,:,:)); % [k x p]
    lam(i,:) = exp(Xi*gamma_pred(:,i) + offset(i,:)')';
    Sigma_pred(:,:,i) = F*Sigma(:,:,i-1)*F' + Q;

    obsidx = isfinite(n(i,:)) & isfinite(lam(i,:));
    
    Wpostinv = inv(Sigma_pred(:,:,i)) + Xi(obsidx,:)'*diag(lam(i,obsidx))*Xi(obsidx,:); % [p x tau][tau x tau][tau x p]
    Sigma(:,:,i) = inv(Wpostinv);

    gamma(:,i)  = gamma_pred(:,i) + Sigma(:,:,i)*Xi(obsidx,:)'*(n(i,obsidx)-lam(i,obsidx))';

    if any(~isfinite(lam(i,:)))
        keyboard
    end
    [~, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:illConditionedMatrix')
        fprintf('MATLAB:illConditionedMatrix at Iter %i \n',i)
        keyboard
        return;
    end
end
toc


predNLL = -sum(n(:).*log(lam(:))-lam(:));

if doSmooth
    fprintf('Smoothing...')
    tic
    I = eye(size(X,3));
    % Backward-Pass (RTS)
    for i=(size(n,1)-2):-1:1
        Sigmai = inv(Sigma_pred(:,:,i+1));
        Fsquig = inv(F)*(I-Q*Sigmai);
        Ksquig = inv(F)*Q*Sigmai;

        gamma(:,i)=Fsquig*gamma(:,i+1) + Ksquig*gamma_pred(:,i+1);
        C = Sigma(:,:,i)*F'*Sigmai;
        Sigma(:,:,i) = Sigma(:,:,i) + C*(Sigma(:,:,i+1)-Sigma_pred(:,:,i+1))*C';
        lam(i,:) = exp(squeeze(X(i,:,:))*gamma(:,i) + offset(i,:)')';
    end
    toc
else
    % recalculate lamda to reflect local likelihoods
    for i=(size(n,1)-2):-1:1
        lam(i,:) = exp(squeeze(X(i,:,:))*gamma(:,i) + offset(i,:)')';
    end
end
