
function mdl = fit_linear_model(data,ij,params,mdl)
% assumes time and Tlist start at t=0 and that time and lfp are regularly sampled


if nargin<4
    mdl.corr_type='none';
    mdl.chan = 1;
    mdl.covariate_fn = @(x) covariate_generator(x);
    mdl.train_id_range = [1 size(data.lfp,1)];
end

% % adjusting spike time to lfp!
% n1 = zeros(size(data.lfp(mdl.train_id_range(1):mdl.train_id_range(2),:),1),1);
% n1(ceil(data.Tlist{ij(1)}(data.Tlist{ij(1)}<data.time(mdl.train_id_range(2)))/params.sta_dt))=1;
% n2 = zeros(size(n1));
% n2(ceil(data.Tlist{ij(2)}(data.Tlist{ij(2)}<data.time(mdl.train_id_range(2)))/params.sta_dt))=1;

n1 = histcounts(data.Tlist{ij(1)},(data.time(mdl.train_id_range(1))-params.sta_dt/2):params.sta_dt:(data.time(mdl.train_id_range(2))+params.sta_dt));
n2 = histcounts(data.Tlist{ij(2)},(data.time(mdl.train_id_range(1))-params.sta_dt/2):params.sta_dt:(data.time(mdl.train_id_range(2))+params.sta_dt));
n1=n1'; n2=n2';

mdl.X = mdl.covariate_fn(data.lfp(mdl.train_id_range(1):mdl.train_id_range(2),mdl.chan));

if nargin<4 || ~isfield(mdl,'beta1')
    % estimate beta
    mdl.beta1 = inv(mdl.X'*mdl.X)*mdl.X'*n1;
    mdl.beta2 = inv(mdl.X'*mdl.X)*mdl.X'*n2;
end


% linear predictor of lfp
mdl.a = xcorr(mdl.X,params.sta_nbins,mdl.corr_type);
mdl.c = xcorr(n2,n1,params.sta_nbins,mdl.corr_type);
mdl.chat=[];

b1 = kron(eye(size(mdl.X,2)),mdl.beta1);
mdl.chat = mdl.a*b1*mdl.beta2;
mdl.c(mdl.c<0)=0;
[mdl.g,mdl.dv] = glmfit(mdl.chat,mdl.c,'poisson');
mdl.gchat = glmval(mdl.g, mdl.chat, 'log');