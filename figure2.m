%figure 2
clc;
clear all;
close all;
%%
load('ecephys715093703_sample');

%%
dt=1./1250;
n1 = zeros(1,length(data.time))';
n1(ceil(data.Tlist{1}/dt))=1;
n2 = zeros(1,length(data.time))';
n2(ceil(data.Tlist{3}/dt))=1;
tid = ceil(529/dt):ceil(531/dt);
chan=1;
spiketrain1 = n1(tid);
spiketrain2 = n2(tid);
%% A) real LFP trace, real spikes from two example neurons
% plot simulated data
figure(1);
% xl=[max(time(ceil(529/dt))) max(time(ceil(531/dt)))];
subplot(3,1,1);
plot(data.time(tid),data.lfp(tid,chan));
% xlim(xl)
grid on;
set(gca,'TickDir','out')
subplot(3,1,2);
stem(data.time(tid), spiketrain1 , 'r', 'Marker', 'none');
% xlim(xl)
set(gca,'TickDir','out')
subplot(3,1,3);
stem(data.time(tid), spiketrain2, 'b', 'Marker', 'none');
% xlim(xl)
xlabel('Time (s)');
set(gca,'TickDir','out')
% exportgraphics(gcf, 'fig2a.pdf', 'ContentType', 'vector');

%% B & C) Frequency model – STAs and power spectra for real data
% fit models
fs=1250;
params.sta_nbins = fs/2;
params.sta_t = linspace(-params.sta_nbins*params.sta_dt,params.sta_nbins*params.sta_dt,2*params.sta_nbins+1);
mdl=[];
mdl.corr_type = 'none';
mdl.covariate_fn = @(x) covariate_generator(x);
mdl.chan =2;%
mdl.maxlag=2000/2;
mdl.Fc = 60;
ij=[1 3];
mdl.train_id_range = [1 800000]; % 1200000 hippo
mdl_freq1 = fit_freq_model(data, ij, params, mdl);

%% C) plot correlograms
figure(2);
clf
bar(params.sta_t,mdl_freq1.c,1);
hold on
plot(params.sta_t,mdl_freq1.gchat)
% ylim([0 0.7]);
hold off
legend('cross correlation','model')
% exportgraphics(gcf, 'fig2b_3.pdf', 'ContentType', 'vector');
%% B)
% plot STAs
figure();
clf
plot(params.sta_t,mdl_freq1.sta1)
hold on
plot(params.sta_t,mdl_freq1.sta2)
axis tight
hold off
legend('STA1', 'STA2')
% ylim([-7e-5 7e-5])
set(gca,'TickDir','out')

% exportgraphics(gcf, 'fig2b_1.pdf', 'ContentType', 'vector');


% plot LFP spectra
f = (0:length(mdl_freq1.a))/length(mdl_freq1.a)*1/params.sta_dt/2;
figure(4);
clf
f = f(f<mdl.Fc);
plot(f,abs(mdl_freq1.Pa));
xlim([0 30]);
set(gca,'TickDir','out')
exportgraphics(gcf, 'fig2b_2.pdf', 'ContentType', 'vector');
%% D & E & F) trace of the generalized phase
fs=1250;
params.sta_nbins = fs/2;
params.sta_t = linspace(-params.sta_nbins*params.sta_dt,params.sta_nbins*params.sta_dt,2*params.sta_nbins+1);
mdl=[];
mdl.chan = 1;
mdl.corr_type = 'unbiased';
mdl.covariate_fn = @(x) covariate_generator(x,'gp');
mdl.cov_mode = 'win';

mdl.train_id_range = [1 800000]; % 1200000 for hippo
mdl_gp = fit_linear_model(data,ij ,params,mdl);
%% D & F) plot correlograms
figure(5);
clf
bar(params.sta_t,mdl_gp.c,1);
hold on
plot(params.sta_t,mdl_gp.gchat)
% ylim([0 0.7]);
hold off
exportgraphics(gcf, 'fig2d_3.pdf', 'ContentType', 'vector');

% plot trace
figure(6);
subplot(3,1,1)
plot(data.time(tid),mdl_gp.X(tid,2));
hold on;
plot(data.time(tid),mdl_gp.X(tid,3),'--');
legend('band pass filtered signal(real)','imag')


c=mdl_gp.X(:,2)+1i*mdl_gp.X(:,3);
cr=ceil((angle(c)+pi)/2/pi*256);
cmap=hsv(256);
subplot(3,1,2)
scatter(data.time(tid),mdl_gp.X(tid,2),20,cmap(cr(tid),:),'filled')
axis tight
% exportgraphics(gcf, 'fig2c.pdf', 'ContentType', 'vector');
%% E) Bilinear model (GP) - betas and cross-covariance of feature
figure(8);
series = [5, 8, 6, 9];
j=0;
subplot(2, 2, 1);
plot(mdl_gp.a(:, series(1)));
title('Signal 2,2');
ylim([-6.5e-9 6.5e-9])
xlim([0 1251])
set(gca,'TickDir','out')

subplot(2, 2, 2);
plot(params.sta_t,mdl_gp.a(:, series(2)));
title('Signal 2,3');
ylim([-6.5e-9 6.5e-9])
set(gca,'TickDir','out')


subplot(2, 2, 3);
plot(params.sta_t,mdl_gp.a(:, series(3)));
title('Signal 3,2');
ylim([-6.5e-9 6.5e-9])
set(gca,'TickDir','out')

subplot(2, 2, 4);
plot(params.sta_t,mdl_gp.a(:, series(4)));
title('Signal 3,3');
ylim([-6.5e-9 6.5e-9])
set(gca,'TickDir','out')
% exportgraphics(gcf, 'fig2d_1.pdf', 'ContentType', 'vector');

% plot Beta
figure(10);
clf
Beta=[mdl_gp.beta1';mdl_gp.beta2']';
h = bar(Beta);
% transparency
h(1).FaceAlpha = 0.8;
h(2).FaceAlpha = 0.8;
legend('Beta 1','Beta2')
% exportgraphics(gcf, 'fig2d_2.pdf', 'ContentType', 'vector');

%% G&I) trace of filterbank outputs (as a heatmap/wavelet transform)
mdl=[];
mdl.chan = 2;
mdl.corr_type = 'unbiased';
mdl.covariate_fn = @(x) covariate_generator(x,'fbank_lowd');
mdl.train_id_range = [1 1000000];
mdl.cov_mode = 'win';

mdl_fb = fit_linear_model(data,ij,params,mdl);
%
% plot correlograms
figure(11);
clf
bar(params.sta_t,mdl_fb.c,1);
hold on
plot(params.sta_t,mdl_fb.gchat)
hold off

% exportgraphics(gcf, 'fig2f_3.pdf', 'ContentType', 'vector');

figure(12);
% plot trace
imagesc(mdl_fb.X(tid,2:9)')
set(gca,'TickDir','out')
colorbar;  
% exportgraphics(gcf, 'fig2E.pdf', 'ContentType', 'vector');
%%  H) Bilinear model (filterbank) - betas and cross-covariance
figure();
%(4,4)  (6,6) (9,9)
series = [52,86,137];

subplot(3, 3, 1);
plot(params.sta_t,mdl_fb.a(:, series(1)));
title('Signal 2,2(0,0)');
ylim([-3*10^-8 3*10^-8])
% xlim([0 1251])
set(gca,'TickDir','out')

subplot(3, 3, 2);
plot(params.sta_t,mdl_fb.a(:, series(2)));
title('Signal 9,2(64,0)');
ylim([-3*10^-8 3*10^-8])
set(gca,'TickDir','out')

subplot(3, 3, 3);
plot(params.sta_t,mdl_fb.a(:, series(3)));
title('Signal 9,9(64,64)');
ylim([-3*10^-8 3*10^-8])
set(gca,'TickDir','out')

% exportgraphics(gcf, 'fig2f_1.pdf', 'ContentType', 'vector');
% plot Beta
figure(14);
clf
plot(mdl_fb.beta1,'-o')
hold on;
plot(mdl_fb.beta2,'-o')
xticks([1:16])
set(gca,'TickDir','out')
legend('Beta 1','Beta2')
% exportgraphics(gcf, 'fig2f_2l.pdf', 'ContentType', 'vector');