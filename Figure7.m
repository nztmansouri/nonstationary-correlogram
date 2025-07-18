
%% dependences
% https://github.com/mullerlab/wave-matlab
% https://github.com/mullerlab/generalized-phase
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
% https://ltfat.org/doc/filterbank/
%%

%% data params

params.data_dir = 'Z:\data\';
params.ecephys_id = '715093703';
params.eceprobe_id = '810755797'; % hc example
% params.eceprobe_id = '810755801'; % ctx example
params.sta_dt = 1/1250;
%
data = load_data_area(params,'hc');


data = load_data_area(params,'hc');

%%  load specific example

load('recording_list.mat')
% load('rec_pair_list_ctx.mat')
% rec=14 ;% 14 &17 (10,26) ;%47 & 19 >>(4,5);
% file_name = ctx_rec_list(rec).name;

load('rec_pair_list_hc.mat')
rec=42 ; % /25,21 (9 13) /42 , 7 (3,15)
file_name = hc_rec_list(rec).name;


% Extract numeric parts
numbers = regexp(file_name, '\d+', 'match');
params.ecephys_id = numbers{1};
params.eceprobe_id = numbers{2}; % hc example
params.sta_dt = 1/1250;

data = load_data_area(params,'hc');
% data = load_data_area(params,'ctx');
params.T = max(cellfun(@max,data.Tlist));


%%

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


%% fit linear models for predicting spike-spike xcorr

ij = [4 5]; % pair from Fig 6
% ij = [5 12];
ij = rec_pair_list{rec}(7,:)

% check pair details
data.units_subset(ij,'ecephys_structure_acronym')
data.units_subset{ij(1),'probe_id'}
[chan_diff,closest_chan]=min(abs(data.chans-data.units_subset{ij(1),'peak_channel_id'}))


mdl=[];
mdl.chan = [closest_chan];
mdl.corr_type = 'none';
mdl.covariate_fn = @(x) covariate_generator(x,'gp');
% mdl.train_id_range = [1 1250*70*60];
mdl.train_id_range = [1250*60*60 1250*70*60];
mdl.cov_mode = 'win';

mdl = fit_linear_model(data,ij,params,mdl);


%% single example

params.window_size=10;
params.window_stride=10;
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



[fsc_no_offset,dsc_no_offset] = glmfit(zscore(cwhat(:)),cwobs(:),'poisson');
% gcwhat0=glmval(fsc,cwhat(:),'log');
% gcwhat0=reshape(gcwhat0,size(cwhat));



figure(16)
ax=[];
ax(1)=subplot(3,1,1);
imagesc(t/60,params.t,cwobs_plot)
cl=get(gca,"CLim");
box off; set(gca,'TickDir','out')
ylabel('Data')
ylim([-.5 .5])
colorbar
ax(2)=subplot(3,1,2);
imagesc(t/60,params.sta_t,gcwhat)
set(gca,'Clim',cl/(size(cwobs,1)/params.tn))
box off; set(gca,'TickDir','out')
ylabel('Nonstationary Model')
colorbar
ax(3)=subplot(3,1,3);
imagesc(t/60,params.sta_t,gcwhat_stationary)
% imagesc(t/60,params.sta_t,gcwhat0)
set(gca,'Clim',cl/(size(cwobs,1)/params.tn))
ylabel('Stationary Model')
xlabel('Time [min]')
colorbar
linkaxes(ax)
xlim([1 45])

[(d0-dsc)/d0 (d0-dsc_stationary)/d0]

%% Basic adaptive smoother demo


b0 = fsc_no_offset;
W0 = eye(length(mdl.g)); % initial parameter estimates
F = eye(length(mdl.g)); % linear dynamics
% Q = eye(length(mdl.g))*10e-2; % process noise
Q = diag([10e-3 10e-8]); % process noise
W0=Q;

cwhat_norm = cwhat-mean(cwhat(:));
cwhat_norm = cwhat_norm./std(cwhat_norm(:));
X = cat(3,ones(size(cwhat')),cwhat_norm'); % add intercept

% [predLL,gamma,Sigma,lam] = ppasmoo_poissexp_vecobs(cwobs',X,b0,W0,F,Q,[],false); % filter only
[predLL,gamma,Sigma,lam] = ppasmoo_poissexp_vecobs(cwobs',X,b0,W0,F,Q,[],true);

cw_asmoo = lam';
dev_asmoo = 2*(sum(cwobs(:).*log(cwobs(:)+(cwobs(:)==0))-cwobs(:))-sum(cwobs(:).*log(cw_asmoo(:))-cw_asmoo(:)));
(d0-dev_asmoo)/d0


%
figure(17)
clf
ax=[];
ax(1)=subplot(5,1,1);
% plot(t/60,gamma(1,:)')
shadedErrorBar(t/60,gamma(1,:)',sqrt(squeeze(Sigma(1,1,:))))
axis tight
ylabel('Baseline')
ax(2)=subplot(5,1,2);
% plot(t/60,gamma(2,:)')
shadedErrorBar(t/60,gamma(2,:)',sqrt(squeeze(Sigma(2,2,:))))
axis tight
ylabel('Gain')
ax(3)=subplot(5,1,3);
imagesc(t/60,params.t,cwobs_plot)
ylabel('Data')
ax(4)=subplot(5,1,4);
imagesc(t/60,params.t,lam')
ylabel('Adaptive Smoother')
clim([min(gcwhat(:)) max(gcwhat(:))])
%
% lam_lfp = exp(squeeze(X(:,:,2)).*gamma(2,:)');
% imagesc(t/60,params.t,lam_lfp')
ax(5)=subplot(5,1,5);
imagesc(t/60,params.t,gcwhat)
ylabel('Indirect Model')
linkaxes(ax,'x')
% xlim([0 45])


%% Optimize Q by grid search on a small set of data

q1=linspace(-9,-2,6)
q2=linspace(-9,-2,6);
for i=1:length(q1)
    for j=1:length(q2)
        QQ = diag(10.^[q1(i) q2(j)]);
        [predNLL,gamma,Sigma,lam,] = ppasmoo_poissexp_vecobs(cwobs(:,1:500)',X(1:500,:,:),b0,W0,F,QQ,[],false);
        predNLLgrid(i,j) = predNLL;
    end
end

%% Optimize Q
fopt = @(b) ppasmoo_poissexp_vecobs(cwobs(:,1:500)',X(1:500,:,:),b0,W0,F,diag(exp(b)),[],false);

options=optimset('Display','iter','MaxFunEvals',100);
bopt = fminsearch(fopt,log(10.^[-3 -6]),options)
Qopt = diag(exp(bopt));

%%


b0 = fsc; W0 = eye(length(mdl.g)); % initial parameter estimates
F = eye(length(mdl.g)); % linear dynamics
% Q = eye(length(mdl.g))*10e-3; % process noise
Q=Qopt;

X = cat(3,ones(size(cwhat')),cwhat'); % add intercept

[predLL,gamma,Sigma,lam] = ppasmoo_poissexp_vecobs(cwobs',X,b0,W0,F,Q,[],true);
cw_asmoo = lam';
dev_asmoo = 2*(sum(cwobs(:).*log(cwobs(:)+(cwobs(:)==0))-cwobs(:))-sum(cwobs(:).*log(cw_asmoo(:))-cw_asmoo(:)));
(d0-dev_asmoo)/d0


%
figure(7)
clf
ax=[];
ax(1)=subplot(5,1,1);
% plot(t/60,gamma(1,:)')
shadedErrorBar(t/60,gamma(1,:)',sqrt(squeeze(Sigma(1,1,:))))
axis tight
ylabel('Baseline')
ax(2)=subplot(5,1,2);
% plot(t/60,gamma(2,:)')
shadedErrorBar(t/60,gamma(2,:)',sqrt(squeeze(Sigma(2,2,:))))
axis tight
ylabel('Gain')
ax(3)=subplot(5,1,3);
imagesc(t/60,params.t,cwobs_plot)
ax(4)=subplot(5,1,4);
imagesc(t/60,params.t,lam')
clim([min(gcwhat(:)) max(gcwhat(:))])
%
% lam_lfp = exp(squeeze(X(:,:,2)).*gamma(2,:)');
% imagesc(t/60,params.t,lam_lfp')


ax(5)=subplot(5,1,5);
imagesc(t/60,params.t,gcwhat)
linkaxes(ax,'x')
xlim([1 5])

%% illustrating flexibility...


b0 = [2.5 0.2]'; W0 = eye(length(mdl.g)); % initial parameter estimates
F = eye(length(mdl.g)); % linear dynamics
% Q = eye(length(mdl.g))*10e-2; % process noise
Q = diag([10e-5 10e-6]); % process noise
W0=Q;

cwhat_norm = cwhat-mean(cwhat(:));
cwhat_norm = cwhat_norm./std(cwhat_norm(:));
cwhat_norm = circshift(cwhat_norm',250)';
X = cat(3,ones(size(cwhat')),cwhat_norm'); % add intercept

% [predLL,gamma,Sigma,lam] = ppasmoo_poissexp_vecobs(cwobs',X,b0,W0,F,Q,[],false); % filter only
% [predLL,gamma,Sigma,lam] = ppasmoo_poissexp_vecobs(cwobs',X,b0,W0,F,Q,[],true);


gamma0_true = 2.5*ones(size(cwhat_norm,2),1);
gamma1_true = 0.2*ones(size(cwhat_norm,2),1);
gamma1_true(150:160) = 0.25;
gamma1_true(161:170) = 0.2;
gamma1_true(171:180) = 0.25;
gamma1_true(181:190) = 0.2;
gamma1_true(191:200) = 0.25;
lam_true = exp(gamma0_true'+gamma1_true'.*cwhat_norm);
% cwobstmp(:,150:200)=poissrnd(exp(1.5+0.2*cwhat_norm(:,150:200)));
cwobstmp = poissrnd(lam_true);
cwobstmp(:,50:100)=NaN; % test missing data

[predLL,gamma,Sigma,lam] = ppasmoo_poissexp_vecobs(cwobstmp',X,b0,W0,F,Q,[],false);

cw_asmoo = lam';
dev_asmoo = 2*(sum(cwobs(:).*log(cwobs(:)+(cwobs(:)==0))-cwobs(:))-sum(cwobs(:).*log(cw_asmoo(:))-cw_asmoo(:)));
(d0-dev_asmoo)/d0


%
figure(7)
clf
ax=[];
ax(1)=subplot(5,1,1);
% plot(t/60,gamma(1,:)')
shadedErrorBar(t/60,gamma(1,:)',sqrt(squeeze(Sigma(1,1,:))))
axis tight
ylabel('Baseline')
hold on
line(xlim(),[1 1]*b0(1))
plot(t/60,gamma0_true)
legend('estimate','initial value','true baseline');
hold off
ax(2)=subplot(5,1,2);
% plot(t/60,gamma(2,:)')
shadedErrorBar(t/60,gamma(2,:)',sqrt(squeeze(Sigma(2,2,:))))
axis tight
ylabel('Gain')
hold on
line(xlim(),[1 1]*b0(2))
plot(t/60,gamma1_true)

hold off
legend('estimate','initial value',' true Gain');


ax(3)=subplot(5,1,3);
% imagesc(t/60,params.t,cwobs_plot)
imagesc(t/60,params.t,cwobstmp)
ylabel('Data')
ax(4)=subplot(5,1,4);
imagesc(t/60,params.t,lam')
ylabel('Adaptive Smoother')
clim([min(gcwhat(:)) max(gcwhat(:))])
clim([min(cwobstmp(:)) max(cwobstmp(:))])
colorbar;                        % Add colorbar
  

% lam_lfp = exp(squeeze(X(:,:,2)).*gamma(2,:)');
% imagesc(t/60,params.t,lam_lfp')
ax(5)=subplot(5,1,5);
% imagesc(t/60,params.t,gcwhat)
% ylabel('Indirect Model')
imagesc(t/60,params.t,lam_true)
ylabel('True Rate')
linkaxes(ax,'x')
xlim([0 45])
%%
%Cortex examples: a) rec_pair_list{47}(19,:)