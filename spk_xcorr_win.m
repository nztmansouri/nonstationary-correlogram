function [corr_win,n,t] = spk_xcorr_win(Tlist,ij,params,corr_type)

% T = max(cellfun(@max,Tlist));
T = params.T;
corr_win = zeros(params.tn,ceil((T-params.window_size)/params.window_stride));
n = zeros(ceil((T-params.window_size)/params.window_stride),2);

% start iteration
window_start = 0;
j=1;
tic
while (window_start+params.window_size) <= T
    fprintf('.')
    window_end = window_start + params.window_size;
    window_Tlist = cellfun(@(x) x(x>window_start & x<window_end),Tlist,'UniformOutput',false);

    switch lower(corr_type)
        case 'both'
            if ~(isempty(window_Tlist{ij(1)}) || isempty(window_Tlist{ij(2)}))
                corr_win(:,j) = spk_xcorr(window_Tlist{ij(1)},window_Tlist{ij(2)},params.trange(1),params.trange(2),params.tn);
            end
        case '1'
            if ~(isempty(window_Tlist{ij(1)}) || isempty(Tlist{ij(2)}))
                corr_win(:,j) = spk_xcorr(window_Tlist{ij(1)},Tlist{ij(2)},params.trange(1),params.trange(2),params.tn);
            end
        case '2'
            if ~(isempty(Tlist{ij(1)}) || isempty(window_Tlist{ij(2)}))
                corr_win(:,j) = spk_xcorr(Tlist{ij(1)},window_Tlist{ij(2)},params.trange(1),params.trange(2),params.tn);
            end
        otherwise
            error('unknown correlation type')
    end

    n(j,:) = cellfun(@length,window_Tlist(ij));
    t(j) = (window_start+window_end)/2;

    % Update window start index for the next iteration
    window_start = window_start + params.window_stride;
    j=j+1;
end
toc