
function mdl = fit_bilinear_model(data,ij,params,mdl)
% assumes time and Tlist start at t=0 and that time and lfp are regularly sampled


if nargin<4
    mdl.corr_type='none';
    mdl.chan = 1;
    mdl.covariate_fn = @(x) covariate_generator(x);
    mdl.train_id_range = [1 size(data.lfp,1)];
end
if ~isfield(mdl,'maxiter')
    mdl.maxiter=25;
end
% 
% tij = data.Tlist(ij);
% tij = cellfun(@(x) x(x>(data.time(mdl.train_id_range(1))-params.sta_dt/2)),tij,'UniformOutput',false);
% tij = cellfun(@(x) x(x<(data.time(mdl.train_id_range(2))+params.sta_dt/2)),tij,'UniformOutput',false);
% mdl.c = spk_xcorr(tij{1},tij{2},min(params.sta_t)-params.sta_dt/2,max(params.sta_t)+params.sta_dt/2,params.sta_nbins*2+1);

n1 = histcounts(data.Tlist{ij(1)},(data.time(mdl.train_id_range(1))-params.sta_dt/2):params.sta_dt:(data.time(mdl.train_id_range(2))+params.sta_dt));
n2 = histcounts(data.Tlist{ij(2)},(data.time(mdl.train_id_range(1))-params.sta_dt/2):params.sta_dt:(data.time(mdl.train_id_range(2))+params.sta_dt));
n1=n1'; n2=n2';
mdl.c = xcorr(n2,n1,params.sta_nbins,mdl.corr_type);
mdl.c(mdl.c<0)=0;

mdl.X = mdl.covariate_fn(data.lfp(mdl.train_id_range(1):mdl.train_id_range(2),mdl.chan));
mdl.a = xcorr(mdl.X,params.sta_nbins,'unbiased');


mdl.alph1 = randn(size(mdl.X,2),1)/10;
% mdl.alph1(1) = log(mean(mdl.c));
mdl.alph2 = randn(size(mdl.X,2),1)/10;
b2 = kron(mdl.alph2,eye(size(mdl.X,2)));
[mdl.alph1,mdl.dv] = glmfit(mdl.a*b2,mdl.c,'poisson','constant','off');
mdl.gchat = glmval(mdl.alph1, mdl.a*b2, 'log','constant','off');

dv=[];
for iter=1:mdl.maxiter
    b1 = kron(eye(size(mdl.X,2)),mdl.alph1);
    % mdl.chat = mdl.a*b1*mdl.alph2;

    [mdl.alph2,mdl.dv] = glmfit(mdl.a*b1,mdl.c,'poisson','constant','off');
    mdl.gchat = glmval(mdl.alph2, mdl.a*b1, 'log','constant','off');
    dv=[dv mdl.dv];

    b2 = kron(mdl.alph2,eye(size(mdl.X,2)));
    [mdl.alph1,mdl.dv] = glmfit(mdl.a*b2,mdl.c,'poisson','constant','off');
    mdl.gchat = glmval(mdl.alph1, mdl.a*b2, 'log','constant','off');
    dv=[dv mdl.dv];
end
mdl.dv_trace=dv;

% % [f,dx]=loss_bilinear([mdl.alph1; mdl.alph2],mdl.a,mdl.c);
% loss = @(b) loss_bilinear(b,mdl.a,mdl.c);
% % [diff,g] = autoGrad([mdl.alph1; mdl.alph2],2,f);
% % [dx g]
% keyboard
% 
% options=[];
% options.Method='sd';
% [b,~,~,mf_out] = minFunc(loss,[mdl.alph1;mdl.alph2],options);
% [fopt,dx,lam]=loss_bilinear(b,mdl.a,mdl.c);
% [diff,g] = autoGrad(b,2,loss);
% 
% 
% function [f,dx,lam]=loss_bilinear(a,XX,y)
% 
% a1 = a(1:length(a)/2);
% a2 = a((length(a)/2+1):end);
% 
% b1 = kron(eye(length(a1)),a1);
% b2 = kron(a2,eye(length(a1)));
% xb = XX*b1*a2;
% lam = exp(xb);
% 
% f = -nansum(nansum(y.*xb - lam));
% lam_err = lam-y;
% lam_err(~isfinite(lam_err))=0;
% 
% da1 = (XX*b2)'*(lam_err);
% da2 = (XX*b1)'*(lam_err);
% dx = [da1(:); da2(:)];
% 
% 
