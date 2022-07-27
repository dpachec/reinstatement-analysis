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



%% GO THROUGH FOLDERS RECURSIVELY
% run in folder with subfolders SISP, SIDR, HC, LC
clear, close all

chan2plot = [1 1 1 1 1 1 1 1 1 1 1];
%chan2plot = [2 2 2 2 2 2 2 2 2 2];

tic
folders = dir(); dirs = find(vertcat(folders.isdir));
folders = folders(dirs);


for foldi = 3:length(folders) %start at 3 cause 1 and 2 are . and ...
    
    direct = folders(foldi);
    cd (direct.name)
    sublist = dir('*_rsa.mat');
    sublist = {sublist.name};
    disp (['measurements -> ' num2str(length(sublist))]);
    
    clear allFiles;
    for subji=1:length(sublist)
        load(sublist{subji});
        allFiles{subji,1} = squeeze(rsaZ (chan2plot(subji),:,:,:)); %when more than 1 electrode

    end

    if exist ('timeBins')
        timeBins1 = timeBins (:, [1, end]);
    end


    cd .. % goes up one directory
    %filename = [sublist{subji}(4:end-12)];
    %if chan2plot(1) == 1 str4name = '_HC'; else srt4name = '_LTC';end 
    filename = [folders(foldi).name];
    eval([filename ' = allFiles;']);
    save (filename, filename, '-v7.3');
    %save (filename, filename);
    
    
end


toc


%% plot group data
cond1 = cell2mat(cellfun(@(x) mean(x), SISP, 'un', 0));

d2p = squeeze(mean(cond1, 'omitnan'));

imagesc(d2p)


%set(gca, 'clim', [-.025 .025])

%clustinfo = bwconncomp(h);





%% 







%%
