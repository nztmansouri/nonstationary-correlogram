
function X = covariate_generator(lfp,type)

if nargin<2
    type='h1d';
end

switch lower(type)
    case 'h1d'
        % broadband hilbert transform on single channel
        lfp_h = hilbert(lfp); % hilbert transformation of lfp
        X = [ones(size(lfp,1),1) real(lfp_h) imag(lfp_h)]; % design matrix
    case 'h2d'
        % broadband hilbert transform on dual channel
        lfp_h1 = hilbert(lfp(:,1)); % hilbert transformation of lfp
        lfp_h2 = hilbert(lfp(:,2)); % hilbert transformation of lfp
        X = [ones(size(lfp,1),1) real(lfp_h1) imag(lfp_h1) real(lfp_h2) imag(lfp_h2)]; % design matrix
    case 'hxd'
        X = [ones(size(lfp,1),1)];
        for chan=1:size(lfp,2)
            h = hilbert(lfp(:,chan)); % hilbert transformation of lfp
            X = [X real(h) imag(h)];
        end
    case 'gp'
        % generalized phase method from davis et al. 2000
        gp_parameters.Fs = 1250; % data sampling rate [Hz]
        gp_parameters.filter_order = 4; gp_parameters.f = [5 40]; % filter parameters
        gp_parameters.lp = 0; % cutoff for negative frequency detection [Hz]
    
        x=[];
        x(1,:,:) = lfp';
        xf = bandpass_filter( x, gp_parameters.f(1), gp_parameters.f(2), gp_parameters.filter_order, gp_parameters.Fs );
        xgp = generalized_phase( xf, gp_parameters.Fs, gp_parameters.lp );
        xgp = shiftdim(xgp,1)';
        X = [ones(size(lfp,1),1) real(xgp) imag(xgp)]; % design matrix

        
        % can add additional types here
    case 'fbank'
        [g,a,fc]=waveletfilters(1250, 'bins', 1250, 1, 250, 6, 'morlet','uniform');
        % gf = filterbankfreqz(g,a,1250*10,'plot');
        Y = ufilterbank(lfp,g,1);
        % keyboard
        X = [ones(size(lfp,1),1)];
        for c=1:size(Y,3)
            X = [X real(Y(:,:,c)) imag(Y(:,2:end,c))];
        end
    case 'fbank_lowd'
        [g,a,fc]=waveletfilters(1250, 'bins', 1250, 1, 64, 1, 'morlet','uniform');
        % gf = filterbankfreqz(g,a,1250*10,'plot');
        Y = ufilterbank(lfp,g,1);
        % keyboard
        X = [ones(size(lfp,1),1)];
        for c=1:size(Y,3)
            X = [X real(Y(:,:,c)) imag(Y(:,2:end,c))];
        end
    otherwise
        error('unknown covariate type')
end