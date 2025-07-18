%figure 5
clc;
clear all;
close all;
%%
load('dt.mat');
load('time_4.mat');
load('Tlist_25.mat');
load('Tlist_26.mat');
load('lfp.mat');
%%
data = struct();
data.lfp=lfp;
ij = [1 2]; % all in the Ca1
data.Tlist{ij(1)}=Tlist_25;
data.Tlist{ij(2)}=Tlist_26;
data.time=time;
%%
params = struct();
params.sta_dt=dt;

% correlogram parameters
params.trange = [-0.5 0.5];
params.tn = 202;
params.t = linspace(params.trange(1),params.trange(2),params.tn);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;
params.T = min(max(data.time),max(cellfun(@max,data.Tlist)));

% moving window parameters
params.window_size = 60; % s
params.window_stride = 10; % s

% lfp parameters
params.sta_nbins = 1250/2;
params.sta_dt = 1./1250;
params.sta_t = linspace(-params.sta_nbins*params.sta_dt,params.sta_nbins*params.sta_dt,2*params.sta_nbins+1);
%%
%% fit linear models for predicting spike-spike xcorr

ij = [1 2];

closest_chan=19;


mdl=[];
mdl.chan = [closest_chan];
mdl.corr_type = 'none';
mdl.covariate_fn = @(x) covariate_generator(x,'gp');
mdl.train_id_range = [1 1250*30];
mdl.cov_mode = 'win';

mdl = fit_linear_model(data,ij,params,mdl);
%%
%% single example

params.window_stride=10; % use overlapping windows to make the plot smoother
[cwobs_plot,n,t] = spk_xcorr_win(data.Tlist,ij,params,'1');
params_fine=params;
params_fine.tn=length(params.sta_t);

% get windowed spike-spike correlogram
[cwobs,n,t] = spk_xcorr_win(data.Tlist,ij,params_fine,'1');
% apply model to same windows to get nonstationary predictions
[cwhat,~,a_win] = mdl_xcorr_win(data,params,mdl);

% add output nonlinearity to match linear predictions to data
cwhat_stationary = repmat(mdl.chat,1,size(cwhat,2));



[fsc_stationary,dsc_stationary] = glmfit(cwhat_stationary(:),cwobs(:),'poisson');
gcwhat_stationary=glmval(fsc_stationary,cwhat_stationary(:),'log');
gcwhat_stationary=reshape(gcwhat_stationary,size(cwhat));


offset=repmat(log(n(:,1).*n(:,2)+1)',size(cwhat,1),1);
% [fsc,dsc] = glmfit(cwhat(:),cwobs(:),'poisson');
% gcwhat=glmval(fsc,cwhat(:),'log');
[fsc,dsc] = glmfit(cwhat(:),cwobs(:),'poisson','offset',offset(:));
gcwhat=glmval(fsc,cwhat(:),'log','offset',offset(:));
gcwhat=reshape(gcwhat,size(cwhat));

[~,d0] = glmfit(cwobs(:)*0+1,cwobs(:),'poisson','constant','off');

[fsc,dsc] = glmfit(cwhat(:),cwobs(:),'poisson');
gcwhat0=glmval(fsc,cwhat(:),'log');
gcwhat0=reshape(gcwhat0,size(cwhat));



figure(6)
ax=[];
ax(1)=subplot(3,1,1)
imagesc(t/60,params.t,cwobs_plot)
cl=get(gca,"CLim");
box off; set(gca,'TickDir','out')
ylabel('Data')
ylim([-.5 .5])
colorbar
ax(2)=subplot(3,1,2)
imagesc(t/60,params.sta_t,gcwhat)
set(gca,'Clim',cl/(size(cwobs,1)/params.tn))
box off; set(gca,'TickDir','out')
ylabel('Nonstationary Model')
ax(3)=subplot(3,1,3)
% imagesc(t/60,params.sta_t,gcwhat_stationary)
imagesc(t/60,params.sta_t,gcwhat0)
set(gca,'Clim',cl/(size(cwobs,1)/params.tn))
ylabel('Stationary Model')
xlabel('Time [min]')
linkaxes(ax)
xlim([1 45])

[(d0-dsc)/d0 (d0-dsc_stationary)/d0]
%% LFP Spectrogram
fs=1250;
par.Fs=fs;
par.tapers = [6 11]; % [TW K];
movingwin = [params.window_size, params.window_stride];
[S_lfp,t_lfp,F_lfp] = mtspecgramc(data.lfp(:,closest_chan),movingwin,par);

%%

timel = [ 19.8, 2.7 22.4 , 37.6]; % Example frequencies for the lines

% Define colors (matching your request)
colors = [225, 0.2, 0.8;  % Deep Blue
          0, 0.6, 0.4;  % Emerald Green
          1, 0.75, 0;   % Golden Yellow
          0.8, 0, 0.8]; % Magenta
figure(10);
ax=[];
ax(2)=subplot(4,1,2)
imagesc(t/60,params.t,cwobs_plot)
xlim([1 45])
hold on
% Add horizontal lines at specified 
for i = 1:length(timel)
ax(2)=xline(timel(i), '--r', 'LineWidth', 1);
hold on
end
hold off
cl=get(gca,"CLim");
box off; set(gca,'TickDir','out')
ylabel('Data')
ylim([-.5 .5])

ax(3)=subplot(4,1,3)
imagesc(t/60,params.sta_t,gcwhat)
set(gca,'Clim',cl/(size(cwobs,1)/params.tn))
hold on
% Add horizontal lines at specified 
for i = 1:length(timel)
ax(2)=xline(timel(i), '--r', 'LineWidth', 1);
hold on
end
hold off
box off; set(gca,'TickDir','out')
ylabel('Nonstationary Model')
ax(4)=subplot(4,1,4)
% imagesc(t/60,params.sta_t,gcwhat_stationary)
imagesc(t/60,params.sta_t,gcwhat0)
set(gca,'Clim',cl/(size(cwobs,1)/params.tn))
hold on
% Add horizontal lines at specified 
for i = 1:length(timel)
ax(2)=xline(timel(i), '--r', 'LineWidth', 1);
hold on
end
hold off
ylabel('Stationary Model')
xlabel('Time [min]')
linkaxes(ax)
xlim([1 45])
ax(1)=subplot(4,1,1)
low_freq = F_lfp>0 & F_lfp<60;
imagesc(t_lfp/60,F_lfp(low_freq),(abs(S_lfp(:,low_freq)')));
hold on
% Add horizontal lines at specified 
ylim([0 30 ])
yl = ylim; % Get the current y-axis limits
y_top = yl(2) + (yl(2) - yl(1)) * 0.05; % Position text slightly above the top
for i = 1:length(timel)
ax(1)=xline(timel(i), '--r', 'LineWidth', 1);
    text(timel(i), y_top, sprintf('%.2f', timel(i)), ...
        'VerticalAlignment', 'bottom', ... % Keeps text above the line
        'HorizontalAlignment', 'center', ...
        'Color', 'r', ...
        'FontSize', 7, ...
        'FontWeight', 'bold', ...
        'Clipping', 'off'); % Ensures text stays visible outside the plot
hold on
end
hold off
xlabel('time [min]')
ylabel('Frequency')
set(gca,'YDir','normal')
xlim([1 45])


%%

rescale = mean(diff(params.sta_t))/mean(diff(params.t));
% 19.8 and 2.72 for freq difference
%%22.45 and 37.6 for amplitude diff
figure(11)
clf
subplot(2,2,1:2)
plot(F_lfp(low_freq),abs(S_lfp(round(22.45*60/10)-2,low_freq)))
hold on
plot(F_lfp(low_freq),abs(S_lfp(round(37.6*60/10)-2,low_freq)))
xlim([0 30]);       
legend('time point 22.45 min',' time point 37.6 min')
hold off
subplot(2,2,3)
bar(params.t,cwobs_plot(1:end-1,round(22.45*60/10)-2),1)
hold on
plot(params.sta_t,gcwhat(:,round(22.45*60/10)-2)/rescale)
[~,i]=max(gcwhat(:,round(22.2*60/10)-2));
[~,j]=min(gcwhat(:,round(22.2*60/10)-2));

hold off
yl=ylim();
% line(-[1 1]*2*1/7,ylim(),'Color','r')
% line([1 1]*params.sta_t(i),ylim(),'Color','r')
% line(params.sta_t(i)+2*[[1:4]' [1:4]']*(params.sta_t(i)-params.sta_t(j)),ylim(),'Color','r')
subplot(2,2,4)
% hold on
bar(params.t,cwobs_plot(1:end-1,round(37.6*60/10)-2),1)
hold on
plot(params.sta_t,gcwhat(:,round(37.6*60/10)-2)/rescale)
legend('observed','model')
hold off
ylim(yl)
% line(-[1 1]*2*1/7,ylim(),'Color','r')
% line([1 1]*params.sta_t(i),ylim(),'Color','r')
% line(params.sta_t(i)+2*[[1:4]' [1:4]']*(params.sta_t(i)-params.sta_t(j)),ylim(),'Color','r')

