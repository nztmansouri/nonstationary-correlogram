function [corr_win,t,a_win] = mdl_xcorr_win(data,params,mdl)


p = size(mdl.beta1,1);
b1 = kron(eye(p),mdl.beta1);

if ~isfield(mdl,'cov_mode') || strcmp(mdl.cov_mode,'full')
    X = mdl.covariate_fn(data.lfp(:,mdl.chan));
end

T = params.T;
corr_win = zeros(params.sta_nbins*2+1,ceil((T-params.window_size)/params.window_stride));
a_win = zeros(params.sta_nbins*2+1,p.^2,ceil((T-params.window_size)/params.window_stride));

% start iteration
window_start = 0;
j=1;
tic
while (window_start+params.window_size) <= T
    fprintf('.')
    window_end = window_start + params.window_size;
    tid = data.time>window_start & data.time<=window_end;


    if isfield(mdl,'cov_mode') && strcmp(mdl.cov_mode,'win')
        y = mdl.covariate_fn(data.lfp(tid,mdl.chan));
    else
        y = X(tid,:);
    end

    switch lower(mdl.corr_type)
        case 'none'
            a_win(:,:,j) = xcorr(y,params.sta_nbins);
            corr_win(:,j) = a_win(:,:,j)*b1*mdl.beta2;
        case 'unbiased'
            a_win(:,:,j) = xcorr(y,params.sta_nbins,'unbiased');
            corr_win(:,j) = a_win(:,:,j)*b1*mdl.beta2;
        otherwise
            error('unknown correlation type')
    end
    t(j) = (window_start+window_end)/2;

    % Update window start index for the next iteration
    window_start = window_start + params.window_stride;
    j=j+1;
end
toc