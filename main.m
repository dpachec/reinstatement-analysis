%% main ERS analysis
% example all2all file has shape  (Trials * 2 (Encoding and retrieval) * Chans (1 in this anlysis) * Freqs * Time 

% rsa_m(all2all, window_size, sliding in, frequencies to include, average in time or not)
results = rsa_m(all2all, 500, 100, [1:44], 0);



%% plot

d2p = squeeze(mean(results));
figure()
imagesc(d2p);





%% 