function [EEG] = ref_elec(EEG, montage)

disp ('>>>>> re-referencing...');

chans2exc = cellstr(EEG.chans2exc); 
chans2exc1 = zeros(size(chans2exc,1), 1);
for i = 1:size(chans2exc1, 1)
    chans2exc1(i) = strmatch(chans2exc(i),{EEG.chanlocs.labels}, 'exact');
end

if strcmp (montage, 'aver');
    disp ('>> AVERAGE REFERENCE...');
    dataRef = EEG.data; dataRef(chans2exc1, :) = []; EEG.chanlocs(chans2exc1, :) = [];
    EEG_average = mean(dataRef, 1);
    EEG.data = dataRef - EEG_average;
end


if strcmp (montage, 'bipo');
    disp ('>> BIPOLAR REFERENCE...');
    count = 1;
    EEG.data(chans2exc1, :) = [];EEG.chanlocs(chans2exc1, :) = [];
    currentLetter = 'X'; nextLetter = 'Y';
    clear chan2useList;
    for i = 1:size(EEG.data, 1)-1
        if  currentLetter == nextLetter
            currentLetter = EEG.chanlocs(i).labels(1);
            dataC1 =  squeeze(EEG.data(strcmpi(EEG.chanlocs(i).labels,{EEG.chanlocs.labels}), :,:));
            nextLetter = EEG.chanlocs(i+1).labels(1);
            dataC2 =  squeeze(EEG.data( strcmpi(EEG.chanlocs(i+1).labels,{EEG.chanlocs.labels}), :,:));
            chan2useList{i} = [EEG.chanlocs(i).labels EEG.chanlocs(i+1).labels];
            data = dataC2 - dataC1;
            EEG.data(i,:) = data;%% extract interval between events in the logs
        end 
    end


disp (' >>>> bipolar reference all electrodes');
chan2useList = chan2useList';

end

if strcmp (montage, 'mono');
    EEG.data(chans2exc1, :) = []; EEG.chanlocs(chans2exc1, :) = [];
    disp (' >>>> no re-reference (monopolar)');
end

