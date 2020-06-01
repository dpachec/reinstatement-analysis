%% PROCESS DATA 
clear
load ALLEEG

for si = 1:length(ALLEEG)
    
EEG = ALLEEG{si};    
subj = EEG.subject;

% re-reference the data 
montage = 'bipo'; % 'aver', 'bipo'
EEG = ref_elec(EEG, montage); % EEG

%epoch data into a list of all trials for selected channels
eLim = [-6 6];
[oneListIds oneListData ] = epoch_rec_data(EEG, eLim);


%tf dec
cfg_dec             = [];
cfg_dec.timeRes     = 0.1; % 0.1 = 100ms ; 0.01 = 10ms, or 'all' (brackets) 
cfg_dec.blne        = [-.5 0];  %in secs
cfg_dec.finalCut    = [-4 4];%in secs

[oneListPow] = decompose_rec (oneListData, cfg_dec);

%remove noisy trials 
epo2exc = EEG.epo2exc;
disp('Excluding trials: ')
disp(oneListIds(epo2exc))
oneListPow(epo2exc,:,:,:) = [];
oneListIds(epo2exc) = [];  oneListData(:,:,epo2exc) = [];


end

