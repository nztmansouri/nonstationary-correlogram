
function mdl = fit_freq_model(data,ij,params,mdl)
% assumes time and Tlist start at t=0 and that time and lfp are regularly sampled


if nargin<4
    mdl.corr_type='none';
    mdl.chan = 1;
    mdl.maxlag=2000/2;
    mdl.Fc=60;
    mdl.train_id_range = [1 size(data.lfp,1)];
end


n1 = histcounts(data.Tlist{ij(1)},(data.time(mdl.train_id_range(1))-params.sta_dt/2):params.sta_dt:(data.time(mdl.train_id_range(2))+params.sta_dt/2));
n2 = histcounts(data.Tlist{ij(2)},(data.time(mdl.train_id_range(1))-params.sta_dt/2):params.sta_dt:(data.time(mdl.train_id_range(2))+params.sta_dt/2));
n1=n1'; n2=n2';


% X = zscore(data.lfp(mdl.train_id_range(1):mdl.train_id_range(2),mdl.chan(1)));
X = data.lfp(mdl.train_id_range(1):mdl.train_id_range(2),mdl.chan(1));
                                                                                                                                                                                               

mdl.sta1 = xcorr(n1,X,params.sta_nbins,'none')/sum(n1);
mdl.sta2 = xcorr(n2,X,params.sta_nbins,'none')/sum(n2);

% mdl.sta1 = xcorr(X,n1,params.sta_nbins,'none')/sum(n1);
% mdl.sta2 = xcorr(X,n2,params.sta_nbins,'none')/sum(n2);

% mdl.an1 = xcorr(X*0+1,n1,params.sta_nbins,'none');
% mdl.sta1 = xcorr(X,n1,params.sta_nbins,'none')./mdl.an1;
% mdl.an2 = xcorr(X*0+1,n2,params.sta_nbins,'none');
% mdl.sta2 = xcorr(X,n2,params.sta_nbins,'none')./mdl.an2;
mdl.a = xcorr(X,params.sta_nbins,'unbiased');


mdl.c = xcorr(n2,n1,params.sta_nbins,mdl.corr_type);

% f = 0:30;
f = (0:length(mdl.a))/length(mdl.a)*1/params.sta_dt/2;
f = f(f<mdl.Fc);
fourier_basis = cos(f'*params.sta_t*2*pi) - 1i*sin(f'*params.sta_t*2*pi);
mdl.F1 = fourier_basis*mdl.sta1;
mdl.F2 = fourier_basis*mdl.sta2;
mdl.Pa = fourier_basis*mdl.a;
mdl.chat = real((mdl.F1.*conj(mdl.F2)./mdl.Pa)'*conj(fourier_basis));
mdl.c(mdl.c<0)=0;
[mdl.g,mdl.dv] = glmfit(mdl.chat,mdl.c,'poisson');
mdl.gchat = glmval(mdl.g, mdl.chat, 'log');
 
