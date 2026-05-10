%% loop over all recordings
load('recording_list.mat')
load('rec_pair_list_hc.mat')
pair_list=cell(0);
nsmetrics_all_poiss=cell(0);
nsmetrics_all_jitter=cell(0);

for rec=2: length(hc_features)

    rec
    file_name = hc_rec_list(rec).name;


    % Extract numeric parts
    numbers = regexp(file_name, '\d+', 'match');
    params.ecephys_id = numbers{1};
    params.eceprobe_id = numbers{2}; % hc example
    params.sta_dt = 1/1250;

    data = load_data_area(params,'hc');
    params.T = max(cellfun(@max,data.Tlist));


    %
    rec


    wsizevec = logspace(0,log10(600),32);
    tic
    for w = 1:length(wsizevec)
        params.window_size = wsizevec(w);
        params.window_stride = wsizevec(w);
        % c=1;
        % for i=1:length(data.Tlist)
        % for j=(1+i):length(data.Tlist)
        for c=1:size(rec_pair_list{rec},1)
            i=rec_pair_list{rec}(c,1);
            j=rec_pair_list{rec}(c,2);
            % original data
            [cwobs,n,t] = spk_xcorr_win(data.Tlist,[i j],params,'1');

            % row and col approx
            c_approx = repmat(mean(cwobs),size(cwobs,1),1);
            r_approx = repmat(mean(cwobs,2),1,size(cwobs,2));
            Ls = sum(cwobs(:).*log(cwobs(:)+(cwobs(:)==0))-cwobs(:));
            Dm = 2*(Ls-sum(cwobs(:).*log(c_approx(:)+(c_approx(:)==0))-c_approx(:)));
            D0 = 2*(Ls-sum(cwobs(:).*log(mean(cwobs(:)))-mean(cwobs(:))));
            R2(2) = (D0-Dm)/D0;
            Dm = 2*(Ls-sum(cwobs(:).*log(r_approx(:)+(r_approx(:)==0))-r_approx(:)));
            R2(1) = (D0-Dm)/D0;


            [A,B,R2(3)]=poisson_decomp(cwobs',1);
            rank1_approx = exp(A*B)';
            [A,B,R2(4)]=poisson_decomp(cwobs',2);
            rank2_approx = exp(A*B)';
            R2

            nsmetrics_all_poiss{rec}(c,:,w) = R2;
            % cwobs_all(:,:,c)=cwobs;


            % jittered data
            Rlist = data.Tlist([i j]);
            Rlist{2} = sort(Rlist{2}+randn(size(Rlist{2})));
            [cwobs,n,t] = spk_xcorr_win(Rlist,[1 2],params,'1');

            % row and col approx
            c_approx = repmat(mean(cwobs),size(cwobs,1),1);
            r_approx = repmat(mean(cwobs,2),1,size(cwobs,2));
            Ls = sum(cwobs(:).*log(cwobs(:)+(cwobs(:)==0))-cwobs(:));
            Dm = 2*(Ls-sum(cwobs(:).*log(c_approx(:)+(c_approx(:)==0))-c_approx(:)));
            D0 = 2*(Ls-sum(cwobs(:).*log(mean(cwobs(:)))-mean(cwobs(:))));
            R2(2) = (D0-Dm)/D0;
            Dm = 2*(Ls-sum(cwobs(:).*log(r_approx(:)+(r_approx(:)==0))-r_approx(:)));
            R2(1) = (D0-Dm)/D0;


            [A,B,R2(3)]=poisson_decomp(cwobs',1);
            rank1_approx = exp(A*B)';
            [A,B,R2(4)]=poisson_decomp(cwobs',2);
            rank2_approx = exp(A*B)';
            R2

            nsmetrics_all_jitter{rec}(c,:,w) = R2;

            pair_list{rec}(c,:) = [i j];
            % c=c+1;
            % end
        end
    end

    save('cc_ns_result_hc.mat','nsmetrics_all_jitter','nsmetrics_all_poiss','pair_list')
end
