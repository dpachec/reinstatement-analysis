function [oneListIds oneListData] = ...
    epoch_rec_data (EEG, eLim)

disp ('>> extracting trials based on rooms'); 

chanNames= EEG.chans2plot;
chans2plot = cellstr(EEG.chans2plot);
chans2plot1 = zeros(size(chans2plot,1), 1);
for i = 1:size(chans2plot, 1)
    chans2plot(i);
    chans2plot1(i,:) = strmatch(chans2plot(i),{EEG.chanlocs.labels}, 'exact');
end
%chans2plot1

clear all_enc_ret all_enc_ret1 EEG_b all_enc_retE all_enc_retR;
countE =1; countR =1; countN = 1;
eLim = eLim;
clear rooms itemIDs all_idsR all_idsE;
for ri = 1:4 % rooms
        for i = 1:length(EEG.event)
            if length(EEG.event(i).type) > 1
                %EEG.event(i).type(1)
                if  EEG.event(i).type(1) == num2str(ri-1) & EEG.event(i).type(end) == 'e' 
                    %disp('hola'); %EEG.event(i).type
                    EEG_b = pop_epoch( EEG, {EEG.event(i).type}, eLim, 'newname', 'Continuous EEG Data epochs', 'epochinfo', 'yes');
                    all_enc_retE{countE, :, :} = EEG_b.data(chans2plot1,:);
                    all_idsE{countE} =  EEG.event(i).type;
                    countE = countE+1;
                end
            end
        end
        rooms{ri, 1} = all_enc_retE;
        itemIDs{ri, 1} = all_idsE'; 
        countE = 1; 
        
        for i = 1:length(EEG.event)
            if (length(EEG.event(i).type) > 1)
                if  EEG.event(i).type(3) == num2str(ri-1) & EEG.event(i).type(end) == 'r' 
                    %EEG.event(i).type
                    EEG_b = pop_epoch( EEG, {EEG.event(i).type}, eLim, 'newname', 'Continuous EEG Data epochs', 'epochinfo', 'yes');
                    all_enc_retR{countR, :, :} = EEG_b.data(chans2plot1,:);
                    all_idsR{countR} =  EEG.event(i).type;
                    countR = countR+1;
                end
            end
        end
        countR = 1;
        rooms{ri, 2} = all_enc_retR;
        itemIDs{ri, 2} = all_idsR'; 
             
        clear all_enc_retE all_enc_retR all_idsR all_idsE; 
end

%%add novel items 
for i = 1:length(EEG.event)
    if  str2double(EEG.event(i).type(1)) == 5
        %EEG.event(i).type
        EEG_b = pop_epoch( EEG, {EEG.event(i).type}, eLim, 'newname', 'Continuous EEG Data epochs', 'epochinfo', 'yes');
        all_enc_retN{countN, :, :} = EEG_b.data(chans2plot1,:);
        all_idsN{countN} =  EEG.event(i).type;
        countN = countN+1;
    end
end

oneListRooms = reshape (rooms, 8, 1);
oneListIds = reshape (itemIDs, 8, 1);

oneListIds = vertcat(oneListIds{:}); 
oneListRooms = vertcat(oneListRooms{:});  

%add novel items to the end of the room list
oneListAll = cat(1, oneListRooms, all_enc_retN);
all_idsN = all_idsN';
oneListIds = cat(1, oneListIds, all_idsN);
size(oneListAll)

for i = 1:length(oneListAll)
   oneListData(:,:,i) =  oneListAll{i};
end



disp ('>> trials extracted'); 


