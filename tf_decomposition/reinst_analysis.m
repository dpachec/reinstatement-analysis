%% LOAD DATA (N = 11)

load ('ALLEEG')

%% PROCESS DATA 

for si = 1:length(ALLEEG)

EEG = ALLEEG{si};    
subj = EEG.subject;

%epoch data into a list of all trials for selected channels
eLim = [-6 6];
[oneListIds oneListData ] = epoch_rec_data(EEG, eLim);


%tf dec
cfg_dec             = [];
cfg_dec.timeRes     = 'all'; % 0.1 = 100ms ; 0.01 = 10ms, or 'all' (brackets) 
cfg_dec.blne        = [-.5 0];  %in secs
cfg_dec.finalCut    = [-4 4];%in secs
[oneListPow] = decompose_rec (oneListData, cfg_dec);

%remove noisy trials 
epo2exc = EEG.epo2exc;
disp('Excluding trials: ')
disp(oneListIds(epo2exc))
oneListPow(epo2exc,:,:,:) = []; 
oneListIds(epo2exc) = [];  

%save
format2save = 'pow'; %pow , pha  or rawT
savefile = 1;
chanNames = EEG.chans2plot; 
contr2save = {'SISP' 'SIDR'}; 

[allContrasts] = ...
    exp_all_cond(oneListIds,oneListPow, chanNames, subj,...
    format2save, savefile, contr2save);

end



%% process all data
clearvars; 

%[files] = load_files('*SI_*all2all.mat');
%create_batches(files, 200, 'DI'); %batch bin, condName
[files] = load_files('*all2all.mat');
clearvars 

create_folders('*all2all.mat');

win_size    = 250; 
step_size   = 50; 
onebyone    = 1;    %1 for local rsa, 0 for global RSA
go_through_f(win_size, step_size, onebyone); % this function will go through all folders and generate RSA files inside them








%%
