%%
function [EEG allR] = ...  %% NOT THE SAME FOR ALL (S2 AND THE REST) 
    ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                            tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);

%this script has two exceptions for subject 1 and 8 (juan and josep)
% import for the rest of the subjects is the same

pickResponse = align2resp;
EEG = pop_loadset(datafile_name, datafile_path);
eventChannel = 'TTL';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
EEG.TTL = TTL;

%%extract experiment window only 
EEG.data = EEG.data(:,tEr *EEG.srate:length(EEG.data));
EEG.times = EEG.times(:,tEr*EEG.srate:length(EEG.times));
EEG.pnts = length(EEG.data);

%%import events from data channel
% it is better to use 500 to identify the events (line down here works ok)
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
    'delevent', 'on' );
for i = 1:length(EEG.event) 
   EEG.event(i).type = 'x'; 
end

FID = fopen(datatraj_path, 'r');
if FID < 0, error('Cannot open file'); end
Data = textscan(FID, '%s', 'Delimiter', '\n');

clear times;clear timesEnter;clear timesOut; clear timesRT; clear timesST;
clear stimType;clear stimTest;
Lines = Data{1};%LongLine = sprintf('%s ', Lines{:});
for i = 1:length(Lines)
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';'); 
end
sizes = cellfun(@length,Data); 


count = 0;
for i = 1:length(Data)
   if sizes(i) ~= 3
      count = count+1; 
   end
   if sizes(i) == 3
      rec(i-count,:) = Data{i,1}; % from trajectory file I take time of collision item ID and room
   end
end 

if strcmp(subj, 's01') 
disp ('subj1');

FID = fopen(datalog_path, 'r');
if FID < 0, error('Cannot open file'); end
Data = textscan(FID, '%s', 'Delimiter', '\n');
clear times;clear timesEnter;clear timesOut; clear timesRT; clear timesST;
clear stimType;clear stimTest;
Lines = Data{1};%LongLine = sprintf('%s ', Lines{:});
for i = 1:length(Lines)
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';|'); 
end
sizes = cellfun(@length,Data); 

str = Data{2}{1};B1 = strsplit(str);B1 = B1 (6:85);B1 = B1';
str = Data{3}{1};B2 = strsplit(str);B2 = B2 (6:85);B2 = B2';
for i = 1:size(Data{2})
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';'); 
end

count = 0;

for i = 1:length(Data)
   if sizes(i) == 2 % extract item information
       strTmp = strsplit(string(Data{i}(2)), ':'); 
       items (i,:) = strTmp;
   end
   if sizes(i) ~= 7
      count = count+1; 
   end
   if sizes(i) == 7
      allR(i-count,:) = Data{i,1};
   end
end 

%take position (frame) from block1 for the specific item
items = items(5:84,:); %remove missing files form the beginning
for i = 1:length(items)
   tmp = strsplit(items(i,1));
   items(i,1) = tmp(1); %remove (oldpos:
   tmp = strsplit(items(i,2),')');
   items(i,2) = tmp(1); %remove ') --> X':
end

%%search by names in allR from itmes seen in block2
%to check what was the old frame position of this item (and save it in column 9)
for i = 1:length(items)
    %allR (find (string(allR(:,2))== items(i,1)), 9) = {items(i,1)};
    allR (find (string(allR(:,2))== items(i,1)), 9) = {items(i,2)};
    %method -> %allR (find (string(allR(:,2))== 'zarajo'), 1) = {3}
end

rec (1:80,1) = cellstr(B1);
rec (81:160,1) = cellstr(B2);

%%get the timings at encoding for all old items
recB1 = rec(1:80,:);
for i = 1:length(recB1)
    tmp = str2double(recB1(i,1));
    allR(find(str2double(string(allR(:,9)))==tmp),10) = {recB1(i,3)};
end
allR(:,8) = rec(81:end,3); %% allR > all items at retrieval (column 8 -> time)


% take room from B1 (can be deduced from the frame)

for i = 1 : length(allR) % 80
    tmpStr = str2double(allR{i,9});
    %disp ('start asigning rooms');
    if tmpStr > - 1 && tmpStr < 20
        %disp ('asigning room 0');
        allR(i,11) =  {'_0'}; %%  (allR column 11 -> room at encoding extracted from frame ID)
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if tmpStr >= 20 && tmpStr <= 39
        %disp ('asigning room 1');
        allR(i,11) =  {'_1'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if tmpStr >= 40 && tmpStr <= 59
        %disp ('asigning room 2');
        allR(i,11) =  {'_2'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if tmpStr >= 60 
        allR(i,11) =  {'_3'};
        %disp ('asigning room 3');
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
end

%%> allR > 1)RetPos;2)itemname;3)retrievalRoom;4)userResp;5)reactionTime;
%6)p_np;7)corr_inc;8)timeTriggers;9)oldPosition;10)OldPosTime;11)OldposRoom
%12) congruent_incongruent room


%take confidence H - L 
for i = 1 : length(allR) % 80
    if strcmp(allR{i,4},'1') | strcmp(allR{i,4},'6') 
        allR{i,13} = 'H';
    elseif strcmp(allR{i,4},'2') | strcmp(allR{i,4},'3') | strcmp(allR{i,4},'4') | strcmp(allR{i,4},'5') 
        allR{i,13} = 'L';
    end
end

%take decision times
p_np = allR (:,6); p_np = str2double(p_np);
DTs = allR (:,5); DTs = str2double(DTs);
DTs_pnP = [p_np , DTs]; DTs_pnP (DTs_pnP (:,1) == 2,:) = [];
DTs = DTs_pnP(:,2);
for i = 1 : length(allR) % 80
    if str2double(allR{i,5}) > median (DTs)
        allR{i,14} = {'H'};
    else
        allR{i,14} = {'L'};
    end
end

%take matrix (cong-inc)
for i = 1 : length(allR) % 80
    allR9 = str2double(cellstr(allR{i,9}));
    if str2double(allR{i,1}) == allR9
        %disp('SP')
        allR{i,15} = '0'; %SP
    else
        allR{i,15} = '1'; %DP
    end
end


d = allR (:,8); d = str2double(d);
e = allR (:,10); %e(cellfun('isempty',e(:,1)), :) = {'-1'};  -> useful maybe for another tiem

%%Add events from rec ROOMS (encoding and retrieval): format = e_r_item_e
% FORMAT =
% encodingRoom_retrievalRoom_sameordifferentMatrixPosition_highorLowconf_ ...
%   name_userResponse(1-6)_correctORincorrect(0-1)_eORr
for i = 1:length(d)
if str2double(allR(i,6)) == 2 %% HERE THE NOVEL ITEMS
   EEG.event(length(EEG.event) + 1) = ... %retrieval
        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('5_', allR(i,3),'_' ,allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
end
if strcmp(allR(i,11), '_0') % if encoding room is 0
    if str2double(allR(i,3)) == 0 % if retrieval room is 0
       %disp ('0 and 0 for enc and ret')
       EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 1 
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_1_', allR(i,15), '_',allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 2 
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 3 
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_3_',allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
end
if strcmp(allR(i,11), '_1')
    if str2double(allR(i,3)) == 0
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 1
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 2
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 3
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_3_',allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
end
    if strcmp(allR(i,11), '_2')
        if str2double(allR(i,3)) == 0
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 1
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 2
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 3
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_3_', allR(i,15), '_', allR(i,13), '_',allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
    end

    if strcmp(allR(i,11), '_3')
        if str2double(allR(i,3)) == 0
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_0_', allR(i,15), '_', allR(i,13), '_',allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 1
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type',strcat('3_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_1_',allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 2
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_2_', allR(i,15), '_', allR(i,13),  '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 3
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
    end
   
end

EEG.event = nestedSortStruct(EEG.event, 'latency');
if findCvalue
    EEG = find_closest_values(EEG); % find closest values 
end
EEG = pick_onset_R(EEG, pickResponse)


%%%%%%%%%%% S8

elseif strcmp(subj, 's08')

    disp ('subj8');

FID = fopen(datalog_path, 'r');
if FID < 0, error('Cannot open file'); end
Data = textscan(FID, '%s', 'Delimiter', '\n');
clear times;clear timesEnter;clear timesOut; clear timesRT; clear timesST;
clear stimType;clear stimTest;
Lines = Data{1};%LongLine = sprintf('%s ', Lines{:});
for i = 1:length(Lines)
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';|'); 
end
sizes = cellfun(@length,Data); 

str = Data{2}{1};B1 = strsplit(str);B1 = B1 (6:85);B1 = B1';
str = Data{3}{1};B2 = strsplit(str);B2 = B2 (6:85);B2 = B2';
for i = 1:size(Data{2})
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';'); 
end

count = 0;

for i = 1:length(Data)
   if sizes(i) == 2 % extract item information
       strTmp = strsplit(string(Data{i}(2)), ':'); 
       items (i,:) = strTmp;
   end
   if sizes(i) ~= 7
      count = count+1; 
   end
   if sizes(i) == 7
      allR(i-count,:) = Data{i,1};
   end
end 

%take position (frame) from block1 for the specific item
items = items(5:84,:); %remove missing files form the beginning
for i = 1:length(items)
   tmp = strsplit(items(i,1));
   items(i,1) = tmp(1); %remove (oldpos:
   tmp = strsplit(items(i,2),')');
   items(i,2) = tmp(1); %remove ') --> X':
end

%%search by names in allR from itmes seen in block2
%to check what was the old frame position of this item (and save it in column 9)
for i = 1:length(items)
    %allR (find (string(allR(:,2))== items(i,1)), 9) = {items(i,1)};
    allR (find (string(allR(:,2))== items(i,1)), 9) = {items(i,2)};
    %method -> %allR (find (string(allR(:,2))== 'zarajo'), 1) = {3}
end

rec (1:80,1) = cellstr(B1);
rec (81:160,1) = cellstr(B2);

%%get the timings at encoding for all old items
recB1 = rec(1:80,:);
for i = 1:length(recB1)
    tmp = str2double(recB1(i,1));
    allR(find(str2double(string(allR(:,9)))==tmp),10) = {recB1(i,3)};
end
allR(:,8) = rec(81:end,3); %% allR > all items at retrieval (column 8 -> time)

% take room from B1 (can be deduced from the frame)
for i = 1 : length(allR) % 80
    strTmp = str2double(string(allR(i,9)));
    if strTmp > - 1 && strTmp < 20
        allR(i,11) =  {'_0'}; %%  (allR column 11 -> room at encoding extracted from frame ID)
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if strTmp >= 20 && strTmp <= 39
        allR(i,11) =  {'_1'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if strTmp >= 40 && strTmp <= 59
        allR(i,11) =  {'_2'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if strTmp >= 60 
        allR(i,11) =  {'_3'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
end

%%> allR > 1)RetPos;2)itemname;3)retrievalRoom;4)userResp;5)reactionTime;
%6)p_np;7)corr_inc;8)timeTriggers;9)oldPosition;10)OldPosTime;11)OldposRoom
%12) congruent_incongruent room


%take confidence H - L 
for i = 1 : length(allR) % 80
    if strcmp(allR{i,4},'1') | strcmp(allR{i,4},'6') 
        allR{i,13} = 'H';
    elseif strcmp(allR{i,4},'2') | strcmp(allR{i,4},'3') | strcmp(allR{i,4},'4') | strcmp(allR{i,4},'5') 
        allR{i,13} = 'L';
    end
end

%take decision times
p_np = allR (:,6); p_np = str2double(p_np);
DTs = allR (:,5); DTs = str2double(DTs);
DTs_pnP = [p_np , DTs]; DTs_pnP (DTs_pnP (:,1) == 2,:) = [];
DTs = DTs_pnP(:,2);
for i = 1 : length(allR) % 80
    if str2double(allR{i,5}) > median (DTs)
        allR{i,14} = {'H'};
    else
        allR{i,14} = {'L'};
    end
end

%take matrix (cong-inc)
for i = 1 : length(allR) % 80
    allR9 = str2double(cellstr(allR{i,9}));
    if str2double(allR{i,1}) == allR9
        %disp('SP')
        allR{i,15} = '0'; %SP
    else
        allR{i,15} = '1'; %DP
    end
end


d = allR (:,8); d = str2double(d);
e = allR (:,10); %e(cellfun('isempty',e(:,1)), :) = {'-1'};  -> useful maybe for another tiem


%%Add events from rec ROOMS (encoding and retrieval): format = e_r_item_e
% FORMAT =
% encodingRoom_retrievalRoom_sameordifferentMatrixPosition_highorLowconf_ ...
%   name_userResponse(1-6)_correctORincorrect(0-1)_eORr
for i = 1:length(d)
if str2double(allR(i,6)) == 2 %% HERE THE NOVEL ITEMS
   EEG.event(length(EEG.event) + 1) = ... %retrieval
        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('5_', allR(i,3),'_' ,allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
end
if strcmp(allR(i,11), '_0') % if encoding room is 0
    if str2double(allR(i,3)) == 0 % if retrieavl room is 0
       %disp ('0 and 0 for enc and ret')
       EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 1 
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_1_', allR(i,15), '_',allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 2 
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 3 
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_3_',allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
end
if strcmp(allR(i,11), '_1')
    if str2double(allR(i,3)) == 0
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 1
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 2
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
    if str2double(allR(i,3)) == 3
        EEG.event(length(EEG.event) + 1) = ... %encoding
           struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_3_',allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
       EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
    end
end
    if strcmp(allR(i,11), '_2')
        if str2double(allR(i,3)) == 0
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 1
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 2
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 3
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_3_', allR(i,15), '_', allR(i,13), '_',allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
    end

    if strcmp(allR(i,11), '_3')
        if str2double(allR(i,3)) == 0
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_0_', allR(i,15), '_', allR(i,13), '_',allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_0_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 1
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type',strcat('3_1_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_1_',allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 2
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_2_', allR(i,15), '_', allR(i,13),  '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_2_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if str2double(allR(i,3)) == 3
            EEG.event(length(EEG.event) + 1) = ... %encoding
               struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
           EEG.event(length(EEG.event) + 1) = ... %retrieval
                struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_3_', allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
    end
   
end


EEG.event = nestedSortStruct(EEG.event, 'latency');
if findCvalue
    EEG = find_closest_values(EEG); % find closest values 
end
EEG = pick_onset_R(EEG, pickResponse)


%%%%%%%% rest of the subjects


    
else
disp ('subject other than 1');
clear allR FID Data d e items;

EEG = pop_loadset(datafile_name, datafile_path);
eventChannel = 'TTL';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
%H = strmatch(H, {EEG.chanlocs.labels}, 'exact');
EEG.TTL = TTL;

%%extract experiment window only 
EEG.data = EEG.data(:,tEr *EEG.srate:length(EEG.data));
EEG.times = EEG.times(:,tEr*EEG.srate:length(EEG.times));
EEG.pnts = length(EEG.data);


%%import events from data channel
% it is better to use 500 to identify the events (line down here works ok)
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
    'delevent', 'on' );
for i = 1:length(EEG.event) 
   EEG.event(i).type = 'x'; 
end

%%load trajectory file and log file 
%clear all;
% Extract data from trajectory file 
FID = fopen(datatraj_path, 'r');
if FID < 0, error('Cannot open file'); end

Data = textscan(FID, '%s', 'Delimiter', '\n');
clear times;clear timesEnter;clear timesOut; clear timesRT; clear timesST;
clear stimType;clear stimTest;
Lines = Data{1};%LongLine = sprintf('%s ', Lines{:});
for i = 1:length(Lines)
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';'); 
end
sizes = cellfun(@length,Data); 

count = 0;
for i = 1:length(Data)
   if sizes(i) ~= 5
      count = count+1; 
   end
   if sizes(i) == 5
      rec(i-count,:) = Data{i,1}; % from trajectory file I take time of collision item ID and room
   end
   if sizes(i) == 10 
      recTrajB1(i,:) = Data{i,1};
   end
   if sizes(i) == 11 
      recTrajB2(i,:) = Data{i,1};
   end
end 

%(cellfun('isempty',e(:,1)), :) = {'-1'};  -> useful maybe for another tiem
recTrajB1(cellfun('isempty',recTrajB1(:,1)), :) = [];  
recTrajB2(cellfun('isempty',recTrajB2(:,1)), :) = [];  


%%load log file from recognition setup 
%clear all;
FID = fopen(datalog_path, 'r');
if FID < 0, error('Cannot open file'); end

Data = textscan(FID, '%s', 'Delimiter', '\n');
clear times;clear timesEnter;clear timesOut; clear timesRT; clear timesST;
clear stimType;clear stimTest;
Lines = Data{1};%LongLine = sprintf('%s ', Lines{:});
for i = 1:length(Lines)
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';|'); 
end
sizes = cellfun(@length,Data); 

str = Data{2}{1};B1 = strsplit(str);B1 = B1 (6:85);B1 = B1';
str = Data{3}{1};B2 = strsplit(str);B2 = B2 (6:85);B2 = B2';
for i = 1:size(Data{2})
   Data(i,:) = textscan(Lines{i}, '%s', 'Delimiter', ';'); 
end

count = 0;
for i = 1:length(Data)
   if sizes(i) == 2 % extract item information
       strTmp = strsplit(string(Data{i}(2)), ':'); 
       items (i,:) = strTmp;
   end
   if sizes(i) ~= 9
      count = count+1; 
   end
   if sizes(i) == 9
      allR(i-count,:) = Data{i,1};
   end
end 

items = items(5:84,:); %remove missing files form the beginning
for i = 1:length(items)
   tmp = strsplit(items(i,1));
   items(i,1) = tmp(1); %remove (oldpos:
   tmp = strsplit(items(i,2),')');
   items(i,2) = tmp(1); %remove ') --> X':
end

% search by names in allR from itmes seen in block2
%to check what was the old frame position of this item (and save it in column 9)
for i = 1:length(items)
    %allR (find (string(allR(:,2))== items(i,1)), 9) = {items(i,1)};
    allR (find (string(allR(:,2))== items(i,1)), 9) = {items(i,2)};
    %method -> %allR (find (string(allR(:,2))== 'zarajo'), 1) = {3}
end

rec (1:80,1) = cellstr(B1);
rec (81:160,1) = cellstr(B2);

%%aqui
%%get the timings at encoding for all old items
recB1 = rec(1:80,:);
for i = 1:length(recB1)
    tmp = str2double(recB1(i,1));
    allR(find(str2double(string(allR(:,9)))==tmp),10) = {recB1(i,4)};
end
allR(:,8) = rec(81:end,4); %% allR > all items at retrieval (column 8 -> time)



% take room from B1 (can be deduced from the frame)
for i = 1 : length(allR) % 80
    strTmp = str2double(string(allR(i,9)));
    if strTmp > - 1 && strTmp < 20
        allR(i,11) =  {'_0'}; %%  (allR column 11 -> room at encoding extracted from frame ID)
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if strTmp >= 20 && strTmp <= 39
        allR(i,11) =  {'_1'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if strTmp >= 40 && strTmp <= 59
        allR(i,11) =  {'_2'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
    if strTmp >= 60 
        allR(i,11) =  {'_3'};
        if str2double(allR{i,11}(2)) == str2double(allR(i, 3))
            allR(i,12) =  {'0'};
        else
            allR(i,12) =  {'1'};
        end
    end
end

%%> allR > 1)RetPos;2)itemname;3)retrievalRoom;4)userResp;5)reactionTime;
%6)p_np;7)corr_inc;8)timeTriggers;9)oldPosition;10)OldPosTime;11)OldposRoom
%12) newRoomPos

%take confidence H - L 
for i = 1 : length(allR) % 80
    if strcmp(allR{i,4},'1') | strcmp(allR{i,4},'6') 
        %disp('high');
        allR{i,13} = 'H';
    elseif strcmp(allR{i,4},'2') | strcmp(allR{i,4},'3') | strcmp(allR{i,4},'4') | strcmp(allR{i,4},'5') 
        %disp('low');
        allR{i,13} = 'L';
    end
end

%take decision times
p_np = allR (:,6); p_np = str2double(p_np);
DTs = allR (:,5); DTs = str2double(DTs);
DTs_pnP = [p_np , DTs]; DTs_pnP (DTs_pnP (:,1) == 2,:) = [];
DTs = DTs_pnP(:,2);
for i = 1 : length(allR) % 80
    if str2double(allR{i,5}) > median (DTs)
        allR{i,14} = {'H'};
    else
        allR{i,14} = {'L'};
    end
end

%take matrix (cong-inc)
for i = 1 : length(allR) % 80
    allR9 = str2double(cellstr(allR{i,9}));
    if str2double(allR{i,1}) == allR9
        %disp('SP')
        allR{i,15} = '0'; %SP
    else
        allR{i,15} = '1'; %DP
    end
end



d = allR (:,8); d = str2double(d);
e = allR (:,10); %e(cellfun('isempty',e(:,1)), :) = {'-1'};  -> useful maybe for another tiem


%%Add events from rec ROOMS (encoding and retrieval): format = e_r_item_e
%%Add events from rec ROOMS (encoding and retrieval): format = e_r_item_e
for i = 1:length(d)
        if str2double(allR(i,6)) == 2 %% HERE THE NOVEL ITEMS
            EEG.event(length(EEG.event) + 1) = ... %retrieval
            struct('latency', offset + d(i)*EEG.srate, 'type', strcat('5_', allR(i,3),'_' ,allR(i,15), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
        end
        if strcmp(allR(i,11), '_0') % if encoding room is 0
            if str2double(allR(i,3)) == 0 
               EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_0_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_0_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
            if str2double(allR(i,3)) == 1 
                EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_1_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_1_', allR(i,6), '_',allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
            if str2double(allR(i,3)) == 2 
                EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_2_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_2_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
            if str2double(allR(i,3)) == 3 
                EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('0_3_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('0_3_',allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
        end
        if strcmp(allR(i,11), '_1')
            if str2double(allR(i,3)) == 0
                EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_0_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_0_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
            if str2double(allR(i,3)) == 1
                EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_1_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_1_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
            if str2double(allR(i,3)) == 2
                EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_2_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_2_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
            if str2double(allR(i,3)) == 3
                EEG.event(length(EEG.event) + 1) = ... %encoding
                   struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('1_3_',allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
               EEG.event(length(EEG.event) + 1) = ... %retrieval
                    struct('latency', offset + d(i)*EEG.srate, 'type', strcat('1_3_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
            end
        end
            if strcmp(allR(i,11), '_2')
                if str2double(allR(i,3)) == 0
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_0_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_0_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
                if str2double(allR(i,3)) == 1
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_1_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_1_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
                if str2double(allR(i,3)) == 2
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_2_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_2_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
                if str2double(allR(i,3)) == 3
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('2_3_', allR(i,6), '_', allR(i,13), '_',allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('2_3_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
            end

            if strcmp(allR(i,11), '_3')
                if str2double(allR(i,3)) == 0
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_0_', allR(i,6), '_', allR(i,13), '_',allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_0_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
                if str2double(allR(i,3)) == 1
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type',strcat('3_1_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_1_',allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
                if str2double(allR(i,3)) == 2
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_2_', allR(i,6), '_', allR(i,13),  '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_2_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
                if str2double(allR(i,3)) == 3
                    EEG.event(length(EEG.event) + 1) = ... %encoding
                       struct('latency', offset + str2double(e{i})*EEG.srate, 'type', strcat('3_3_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_e'), 'urevent', 0);        
                   EEG.event(length(EEG.event) + 1) = ... %retrieval
                        struct('latency', offset + d(i)*EEG.srate, 'type', strcat('3_3_', allR(i,6), '_', allR(i,13), '_', allR(i,2), '_', allR(i,4), '_', allR(i,7),  '_r'), 'urevent', 0);        
                end
            end
   
end

%%%% include trajectory info as events in EEGLAB
% % if plotTrajEvents
% % %     
% % %     for i = 1:length(recTrajB1)  
% % %         EEG.event(length(EEG.event) + 1) = ... %retrieval
% % %             struct('latency', round (offset + str2double(recTrajB1(i,9))*EEG.srate), 'type', ...
% % %             strcat('t_', recTrajB1(i,1), '_', recTrajB1(i,3), '_', recTrajB1(i,4), '_', ...
% % %             recTrajB1(i,5), '_', recTrajB1(i,6), '_', recTrajB1(i,7), '_', ...
% % %             recTrajB1(i,8)), 'urevent', 0);           
% % %     end
% % %     for i = 1:length(recTrajB2)  
% % %         EEG.event(length(EEG.event) + 1) = ... %retrieval
% % %             struct('latency', round (offset + str2double(recTrajB2(i,10))*EEG.srate), 'type', ...
% % %             strcat('t_', recTrajB2(i,1), '_', recTrajB2(i,4), '_', recTrajB2(i,5), '_', ...
% % %             recTrajB2(i,6), '_', recTrajB2(i,7), '_', recTrajB2(i,8), '_', ...
% % %             recTrajB2(i,9)), 'urevent', 0);      
% % %     end
% % end

%%%% create separate matrix for trajectory info
if plotTrajEvents
    recTrajB1(1:10,:)
    trajInfo = zeros (2, length(EEG.data));
    for i = 1:length(recTrajB1)  
        idx = round (offset + str2double(recTrajB1(i,9))*EEG.srate);
        trajInfo(1, idx) = str2double(recTrajB1(i,3));
        trajInfo(2, idx) = str2double(recTrajB1(i,5));
    end
    for i = 1:length(recTrajB2)  
        idx = round (offset + str2double(recTrajB2(i,10))*EEG.srate);
        trajInfo(1, idx) = str2double(recTrajB2(i,4));
        trajInfo(2, idx) = str2double(recTrajB2(i,6));
    end
    
    EEG.trajInfo = trajInfo;
end

EEG.event = nestedSortStruct(EEG.event, 'latency');


if findCvalue
    EEG = find_closest_values(EEG); % find closest values 
end
EEG = pick_onset_R(EEG, pickResponse);

end

%weird exception
% only for subject 6
if strcmp(subj, 's06')
    EEG.event(2).latency = EEG.event(3).latency;
    
end

disp ('>> Data imported');


%%