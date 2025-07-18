clc;
clear all;
close all;
%%
fs = 1000;      % Sampling frequency (Hz)
dt=1/fs; 

% correlograms and stas...
ntrials=20;
T = 50;
rng(13)

% spike/lfp traces
% ntrials=1;
% T = 5;

data.time = linspace(0, T, T * fs)';  % Time vector


% Simulated LFP #1 - Ramped sin
% f = 4;         % Frequency of the sinusoid (Hz)
% r_rate = 0.25;   % Ramp rate for amplitude increase
% data.lfp = cos(2 * pi * f * data.time) .* ((r_rate * data.time/T*5 + 1));

% % Simulated LFP #2 - Chirp
f0=1; f1=10;
instantaneous_freq = exp((log(f1) - log(f0)) / T * data.time) * f0;
data.lfp = cos(2 * pi * cumsum(instantaneous_freq) / fs);

data.lfp=zscore(data.lfp);

% spike rates
lambda1=1./(1+exp(.75*data.lfp+5.5))*ntrials;
lambda2=1./(1+exp(-.75*data.lfp+5.5))*ntrials;

mean(lambda1)/dt/ntrials

% demo spike trains for two neurons
spike_train1 = poissrnd(lambda1);
spike_train2 = poissrnd(lambda2);


% time of spikes
Tlist=cell(0);

Tlist{1} = data.time(count2list(spike_train1))+dt/2;
Tlist{2} = data.time(count2list(spike_train2))+dt/2;
% jitter to prevent aliasing in correlograms
data.Tlist = cellfun(@(x) x+(rand(size(x))-.5)*dt,Tlist,'UniformOutput',false);



% plot simulated data
figure(1);
xl=[0 max(data.time)];
subplot(3,1,1);
plot(data.time, data.lfp);
xlim(xl)
ylim([-3 3])
grid on;
set(gca,'TickDir','out')
subplot(3,1,2);
stem(data.time, spike_train1, 'r', 'Marker', 'none');
% stem(data.time, spike_train1>=1, 'r', 'Marker', 'none');
xlim(xl)
set(gca,'TickDir','out')
subplot(3,1,3);
stem(data.time, spike_train2, 'b', 'Marker', 'none');
% stem(data.time, spike_train2>=1, 'b', 'Marker', 'none');
xlim(xl)
xlabel('Time (s)');
set(gca,'TickDir','out')


% fit models
ij = [1 2];
params.sta_dt = dt;
params.sta_nbins = fs/2;
params.sta_t = linspace(-params.sta_nbins*params.sta_dt,params.sta_nbins*params.sta_dt,2*params.sta_nbins+1);
mdl=[];
mdl.corr_type = 'none';
mdl.covariate_fn = @(x) covariate_generator(x);
mdl.chan = 1;
mdl.maxlag=2000/2;
mdl.Fc = 60;

mdl.train_id_range =[1 ceil(size(data.lfp,1)/2)]; % first half
mdl_freq1 = fit_freq_model(data, ij, params, mdl);
mdl.train_id_range =[ceil(size(data.lfp,1)/2) size(data.lfp,1)]; % second half
mdl_freq2 = fit_freq_model(data, ij, params, mdl);

% all
mdl.train_id_range =[1 size(data.lfp,1)]; % all
mdl_freq = fit_freq_model(data, ij, params, mdl);


% plot correlograms
figure(2);
clf
ax=[];
ax(1)=subplot(2,2,2);
bar(params.sta_t,mdl_freq2.c,1);
hold on
plot(params.sta_t,mdl_freq2.gchat)
ylim([0 900])
yticks([0 200 400 600 900])
hold off; box off; set(gca,'TickDir','out')
ax(2)=subplot(2,2,1);
bar(params.sta_t,mdl_freq1.c,1);
hold on
plot(params.sta_t,mdl_freq1.gchat)
ylim([0 900])
yticks([0 200 400 600 900])
hold off; box off; set(gca,'TickDir','out')
linkaxes(ax)
legend('spike counts','Model')

% plot STAs

ax=[];
ax(1)=subplot(2,2,4);
plot(params.sta_t,mdl_freq2.sta1)
hold on
plot(params.sta_t,mdl_freq2.sta2)
hold off; box off; set(gca,'TickDir','out')
ax(2)=subplot(2,2,3);
plot(params.sta_t,mdl_freq1.sta1)
hold on
plot(params.sta_t,mdl_freq1.sta2)
hold off; box off; set(gca,'TickDir','out')
linkaxes(ax)
legend('STA1','STA2')




% plot LFP spectra
f = (0:length(mdl_freq1.a))/length(mdl_freq1.a)*1/params.sta_dt/2;
figure(4);
clf
f = f(f<mdl.Fc);
plot(f,abs(mdl_freq1.Pa), f,abs(mdl_freq2.Pa));
ylim([0 1000]);
xlim([0 30]);
hold off; box off; set(gca,'TickDir','out')

figure(3);
bar(params.sta_t,mdl_freq.c,1);
hold on
plot(params.sta_t,mdl_freq.gchat)
ylim([0 900])
yticks([0 200 400 600 800 900])
hold off; box off; set(gca,'TickDir','out')