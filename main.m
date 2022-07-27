
%% main ERS analysis
%% all2all has shape  (Trials * 2 (Encoding and retrieval) * Chans * Freqs * Time 

%all2all= all2all(:, :, 1, :, :);
results = rsa_m(all2all, 500, 100, [1:44], 0);




%% 

d2p = squeeze(mean(results));
figure()
imagesc(d2p);





%% 