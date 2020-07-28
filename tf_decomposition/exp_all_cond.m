

function [allContrasts] = exp_all_cond (oneListIds, eegpower_BN, chanNames, subj, ...
                                        format2save, savefile, contr2save)

allContrasts = 1;

countSISR = 1; countDISR = 1; countSIDR = 1; countDIDR = 1; countSISP = 1; countSIDP = 1;
countSIDPC = 1; countSI = 1; countDIN = 1; countDINC = 1;countDINHCC = 1;
countSR = 1; countDR =1 ; countHC = 1; countLC = 1;
count_HC_Hits = 1; count_HC_Misses = 1; count_HCCRej = 1; count_SIC = 1; count_HCC = 1;
count_SISPC = 1; count_SIDRC = 1; count_DIC = 1;count_SII = 1; count_LCC = 1;
count_1 = 1; count_2 = 1; count_3 = 1; count_4 = 1; count_5 = 1; count_6 = 1; 
count_SISRDINC = 1; count_SIDRDINC = 1; 


%here the main analysis 
for i = 1:length(oneListIds)
    if oneListIds{i}(end) == 'e'
        for j = i:length(oneListIds)
            if oneListIds{j}(end) == 'r' && ... % SI
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) %same item 
               disp (['SI > ' oneListIds{i} '//' oneListIds{j}]);     
               new_si{countSI} = [i, j];
               countSI = countSI+1;
            end
            %here before DI but was not correct cause items were
            %duplicated
            if oneListIds{j}(end) == 'r' && ... % DIN 
               ~strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... % diff item 
               strcmp(oneListIds{j}(1), '5')                           % not presented
               %disp (['DI > ' oneListIds{i} '//' oneListIds{j}]);     
               new_din{countDIN} = [i, j];
               countDIN = countDIN+1;
            end
            if oneListIds{j}(end) == 'r' && ... % DINC %SURROGATES ARE HIGH CONF CRs
               ~strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... % diff item 
               strcmp(oneListIds{j}(1), '5')  && ...                          % not presented
               strcmp(oneListIds{j}(end-2) , '1')                             % correct
               %disp (['DI > ' oneListIds{i} '//' oneListIds{j}]);     
               new_dinc{countDINC} = [i, j];
               countDINC = countDINC+1;
            end
            if oneListIds{j}(end) == 'r' && ... % DINHCC %SURROGATES ARE HIGH CONF CRs
               ~strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... % diff item 
               strcmp(oneListIds{j}(1), '5')  && ...                          % not presented
               strcmp(oneListIds{j}(7), 'H')  && ...                          % HC
               strcmp(oneListIds{j}(end-2) , '1')                             % correct
               %disp (['DI > ' oneListIds{i} '//' oneListIds{j}]);     
               new_dinhcc{countDINHCC} = [i, j];
               countDINHCC = countDINHCC+1;
            end
            if oneListIds{j}(end) == 'r' && ... % SR
               strcmp(oneListIds{i}(1) , oneListIds{j}(3)) && ...  %same room     
               ~strcmp(oneListIds{j}(1) , '5')  %not novel item
               %disp (['SR > ' oneListIds{i} '//' oneListIds{j}]);
               new_sr{countSR} = [i, j];
               countSR = countSR+1;
            end
            if oneListIds{j}(end) == 'r' && ... % DR
               ~strcmp(oneListIds{i}(1) , oneListIds{j}(3)) && ...  %diff room     
               ~strcmp(oneListIds{j}(1) , '5')  %not novel item
               %disp (['DR > ' oneListIds{i} '//' oneListIds{j}]);
               new_dr{countDR} = [i, j];
               countDR = countDR+1;
            end
            if oneListIds{j}(end) == 'r' && ... % SISP
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{j}(5), '0')  %same position
               %disp (['SISP > ' oneListIds{i} '//' oneListIds{j}]);     
               new_sisp{countSISP} = [i, j];
               countSISP = countSISP+1;
            end
            if oneListIds{j}(end) == 'r' && ... % SIDP
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{j}(5), '1')  %diff position
               %disp (['SIDP > ' oneListIds{i} '//' oneListIds{j}]);     
               new_sidp{countSIDP} = [i, j];
               countSIDP = countSIDP+1;
            end          
            if oneListIds{j}(end) == 'r' && ... % SIDPC
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{j}(5), '1')  && ...            %diff position
               strcmp(oneListIds{i}(end-2) , '1') % response is correct
               %disp (['SIDPC > ' oneListIds{i} '//' oneListIds{j}]);     
               new_sidpc{countSIDPC} = [i, j];
               countSIDPC = countSIDPC+1;
            end   
            if oneListIds{j}(end) == 'r' && ... % SISR
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{i}(1) , oneListIds{j}(3))  %same room
               %disp (['SISR > ' oneListIds{i} '//' oneListIds{j}]);     
               new_sisr{countSISR} = [i, j];
               countSISR = countSISR+1;
            end
            if oneListIds{j}(end) == 'r' && ... % DISR
               ~strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... % diff item
               strcmp(oneListIds{i}(1) , oneListIds{j}(3))  %same room
               %disp (['DISR > ' oneListIds{i} '//' oneListIds{j}]);     
               new_disr{countDISR} = [i, j];
               countDISR = countDISR+1;
            end
            if oneListIds{j}(end) == 'r' && ... %SIDR
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                ~strcmp(oneListIds{i}(1) , oneListIds{j}(3)) % diff room
                %disp (['SIDR > ' oneListIds{i} '//' oneListIds{j}]);     
                new_sidr{countSIDR} = [i, j];
                countSIDR = countSIDR+1;
            end
            if oneListIds{j}(end) == 'r' && ... %DIDR
                ~strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... % diff item
                ~strcmp(oneListIds{i}(1) , oneListIds{j}(3)) %% diff room
                %disp (['DIDR > ' oneListIds{i} '//' oneListIds{j}]);     
                new_didr{countDIDR} = [i, j];
                countDIDR = countDIDR+1;
            end
            if oneListIds{j}(end) == 'r' && ... %HC
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                strcmp(oneListIds{i}(7), 'H') %% HC
                %disp (['HC > ' oneListIds{i} '//' oneListIds{j}]);     
                new_hc{countHC} = [i, j];
                countHC = countHC+1;
            end
            if oneListIds{j}(end) == 'r' && ... %LC
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                strcmp(oneListIds{i}(7) , 'L') %% LC
                %disp (['LC > ' oneListIds{i} '//' oneListIds{j}]);     
                new_lc{countLC} = [i, j];
                countLC = countLC+1;
            end
            if oneListIds{j}(end) == 'r' && ... %HC_Hits
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                strcmp(oneListIds{i}(end-4) , '1') && ... %user response is 1
                strcmp(oneListIds{i}(end-2) , '1') % response is correct (dont need this here)
                %disp (['HCC > ' oneListIds{i} '//' oneListIds{j}]);     
                new_hcc{count_HC_Hits} = [i, j];
                count_HC_Hits = count_HC_Hits+1;
            end
            if oneListIds{j}(end) == 'r' && ... %HC_Misses
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                strcmp(oneListIds{i}(end-4) , '6') %user response is 6 (miss)
                %disp (['HC_Hits > ' oneListIds{i} '//' oneListIds{j}]);     
                new_HCI{count_HC_Misses} = [i, j];
                count_HC_Misses = count_HC_Misses+1;
            end
            if oneListIds{j}(end) == 'r' && ... %LC
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                strcmp(oneListIds{i}(7) , 'L') && ... %% LC
                strcmp(oneListIds{i}(end-2) , '1') % response is correct 
                %disp (['LCC > ' oneListIds{i} '//' oneListIds{j}]);     
                new_lcc{count_LCC} = [i, j];
                count_LCC = count_LCC+1;
            end
            if oneListIds{j}(end) == 'r' && ... %SIC
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                strcmp(oneListIds{i}(end-2) , '1') % response is correct
                %disp (['SIC > ' oneListIds{i} '//' oneListIds{j}]);     
                new_SIC{count_SIC} = [i, j];
                count_SIC = count_SIC +1;
            end
            if oneListIds{j}(end) == 'r' && ... %DIC
                ~strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %diff item 
                strcmp(oneListIds{i}(end-2) , '1') && ... % encoding item is correct
                strcmp(oneListIds{j}(end-2) , '1')        % retrieval item is correct
                %disp (['DIC > ' oneListIds{i} '//' oneListIds{j}]);     
                new_DIC{count_DIC} = [i, j];
                count_DIC = count_DIC +1;
            end
            if oneListIds{j}(end) == 'r' && ... %SII
                strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
                strcmp(oneListIds{i}(end-2) , '0') % response is incorrect
                %disp (['SII > ' oneListIds{i} '//' oneListIds{j}]);     
                new_SII{count_SII} = [i, j];
                count_SII = count_SII +1;
            end
            if oneListIds{j}(end) == 'r' && ... % SISPC
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ...  %same item 
               strcmp(oneListIds{j}(5), '0') && ...                           %same position
               strcmp(oneListIds{i}(end-2) , '1')                             %response is correct
               %disp (['SISPC > ' oneListIds{i} '//' oneListIds{j}]);     
               new_SISPC{count_SISPC} = [i, j];
               count_SISPC = count_SISPC+1;
            end
            if oneListIds{j}(end) == 'r' && ... % SIDRC
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               ~strcmp(oneListIds{i}(1) , oneListIds{j}(3)) && ...           % diff room
               strcmp(oneListIds{i}(end-2) , '1')                            %response is correct
               %disp (['SIDRC > ' oneListIds{i} '//' oneListIds{j}]);     
               new_SIDRC{count_SIDRC} = [i, j];
               count_SIDRC = count_SIDRC+1;
            end 
            if oneListIds{j}(end) == 'r' && ... % 1
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{i}(end-4) , '1') %user response is 6 (miss)
               new_1{count_1} = [i, j];
               count_1 = count_1 + 1;
            end 
            if oneListIds{j}(end) == 'r' && ... % 2
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{i}(end-4) , '2') %user response is 2 (miss)
               new_2{count_2} = [i, j];
               count_2 = count_2 + 1;
            end 
            if oneListIds{j}(end) == 'r' && ... % 3
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{i}(end-4) , '3') %user response is 3 (miss)
               new_3{count_3} = [i, j];
               count_3 = count_3 + 1;
            end 
            if oneListIds{j}(end) == 'r' && ... % 4
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{i}(end-4) , '4') %user response is 4 (miss)
               new_4{count_4} = [i, j];
               count_4 = count_4 + 1;
            end 
            if oneListIds{j}(end) == 'r' && ... % 5
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{i}(end-4) , '5') %user response is 5 (miss)
               new_5{count_5} = [i, j];
               count_5 = count_5 + 1;
            end 
            if oneListIds{j}(end) == 'r' && ... % 6
               strcmp(oneListIds{i}(4:end-1), oneListIds{j}(4:end-1)) && ... %same item 
               strcmp(oneListIds{i}(end-4) , '6') %user response is 6 (miss)
               new_6{count_6} = [i, j];
               count_6 = count_6 + 1;
            end 
        end
    end
end

%this is just for the last contrast with novel items as asked for
%the reviewer, has to be done room by room 

%SISRDINC
for i = 1:length(oneListIds)
    if strcmp(oneListIds{i}(1), '0') && strcmp(oneListIds{i}(end-2) , '1') %room 0 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               strcmp(oneListIds{j}(3), '0')        && ...  presented in room 0 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SISRDINC > ' oneListIds{i} '//' oneListIds{j}]); 
               new_SISRDINC{count_SISRDINC} = [i, j];
               count_SISRDINC = count_SISRDINC+ 1;  
            end
        end    
    end
    if strcmp(oneListIds{i}(1), '1') && strcmp(oneListIds{i}(end-2) , '1') %room 1 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               strcmp(oneListIds{j}(3), '1')        && ...  presented in room 1 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SISRDINC > ' oneListIds{i} '//' oneListIds{j}]); 
               new_SISRDINC{count_SISRDINC} = [i, j];
               count_SISRDINC = count_SISRDINC+ 1;  
            end
        end    
    end
    if strcmp(oneListIds{i}(1), '2') && strcmp(oneListIds{i}(end-2) , '1') %room 2 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               strcmp(oneListIds{j}(3), '2')        && ...  presented in room 2 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SISRDINC > ' oneListIds{i} '//' oneListIds{j}]); 
               new_SISRDINC{count_SISRDINC} = [i, j];
               count_SISRDINC = count_SISRDINC+ 1;  
            end
        end    
    end
    if strcmp(oneListIds{i}(1), '3') && strcmp(oneListIds{i}(end-2) , '1') %room 3 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               strcmp(oneListIds{j}(3), '3')        && ...  presented in room 3 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SISRDINC > ' oneListIds{i} '//' oneListIds{j}]); 
               new_SISRDINC{count_SISRDINC} = [i, j];
               count_SISRDINC = count_SISRDINC+ 1;  
            end
        end    
    end
end

%SIDRDINC 
for i = 1:length(oneListIds)
    if strcmp(oneListIds{i}(1), '0') && strcmp(oneListIds{i}(end-2) , '1') %room 0 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               (strcmp(oneListIds{j}(3), '1') || strcmp(oneListIds{j}(3), '2') || ...
               strcmp(oneListIds{j}(3), '3'))      && ...  presented in room 1-2-3 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SIDRDINC > ' oneListIds{i} '//' oneListIds{j}]); 
               new_SIDRDINC{count_SIDRDINC} = [i, j];
               count_SIDRDINC = count_SIDRDINC+ 1;  
            end
        end    
    end
    if strcmp(oneListIds{i}(1), '1') && strcmp(oneListIds{i}(end-2) , '1') %room 1 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               (strcmp(oneListIds{j}(3), '0') || strcmp(oneListIds{j}(3), '2') || ...
               strcmp(oneListIds{j}(3), '3') )     && ...  presented in room 0-2-3 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SIDRDINC > ' oneListIds{i} '//' oneListIds{j}]); 
               new_SIDRDINC{count_SIDRDINC} = [i, j];
               count_SIDRDINC = count_SIDRDINC+ 1;  
            end
        end    
    end
    if strcmp(oneListIds{i}(1), '2') && strcmp(oneListIds{i}(end-2) , '1') %room 2 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               (strcmp(oneListIds{j}(3), '0') || strcmp(oneListIds{j}(3), '1') || ...
               strcmp(oneListIds{j}(3), '3') )     && ...  presented in room 0-1-3 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SIDRDINC > ' oneListIds{i} '//' oneListIds{j}]); 
               new_SIDRDINC{count_SIDRDINC} = [i, j];
               count_SIDRDINC = count_SIDRDINC+ 1;  
            end
        end    
    end
    if strcmp(oneListIds{i}(1), '3') && strcmp(oneListIds{i}(end-2) , '1') %room 3 and correct
        for j = 1:length(oneListIds)
            if strcmp(oneListIds{j}(1), '5')        && ...  novel item
               (strcmp(oneListIds{j}(3), '0') || strcmp(oneListIds{j}(3), '1') || ...
               strcmp(oneListIds{j}(3), '2') )     && ...  presented in room 0-1-2 at retrieval
               strcmp(oneListIds{j}(end-2) , '1')           %is correct
               %disp (['SIDRDINC> ' oneListIds{i} '//' oneListIds{j}]); 
               new_SIDRDINC{count_SIDRDINC} = [i, j];
               count_SIDRDINC = count_SIDRDINC+ 1;  
            end
        end    
    end
end

new_si = new_si'; new_si = vertcat(new_si{:});
new_din = new_din'; new_din = vertcat(new_din{:});
new_dinc = new_dinc'; new_dinc = vertcat(new_dinc{:});
new_dinhcc = new_dinhcc'; new_dinhcc = vertcat(new_dinhcc{:});
new_sr = new_sr'; new_sr = vertcat(new_sr{:});
new_dr = new_dr'; new_dr = vertcat(new_dr{:});
new_sisr = new_sisr'; new_sisr = vertcat(new_sisr{:});
new_sisp = new_sisp'; new_sisp = vertcat(new_sisp{:});
if exist('new_sidp')
new_sidp = new_sidp'; new_sidp = vertcat(new_sidp{:});
else
    new_sidp = {};
end
if exist('new_sidpc')
    new_sidpc = new_sidpc'; new_sidpc = vertcat(new_sidpc{:});
else
   new_sidpc = {};
end
new_disr = new_disr'; new_disr = vertcat(new_disr{:});
if exist('new_sidr')
    new_sidr = new_sidr'; new_sidr = vertcat(new_sidr{:});
else
   new_sidr = {};
end
new_didr = new_didr'; new_didr = vertcat(new_didr{:});
new_hc = new_hc'; new_hc = vertcat(new_hc{:});
if exist('new_lc')
new_lc = new_lc'; new_lc = vertcat(new_lc{:});
else
    new_lc = {};
end
if exist('new_lcc')
new_lcc = new_lcc'; new_lcc = vertcat(new_lcc{:});
else
    new_lcc = {};
end
new_SIC = new_SIC'; new_SIC = vertcat(new_SIC{:});
if exist('new_HCC')
    new_HCC = new_HCC'; new_HCC = vertcat(new_HCC{:});
else
   new_HCC = {}; 
end
new_hcc = new_hcc'; new_hcc = vertcat(new_hcc{:});
if exist('new_HCI')
    new_HCI = new_HCI'; new_HCI = vertcat(new_HCI{:});
else
    new_HCI = {};
end
new_SISPC = new_SISPC'; new_SISPC = vertcat(new_SISPC{:});
if exist('new_SIDRC')
    new_SIDRC = new_SIDRC'; new_SIDRC = vertcat(new_SIDRC{:});
else 
    new_SIDRC ={};
end 
new_DIC = new_DIC'; new_DIC = vertcat(new_DIC{:});
if exist('new_SII')
    new_SII = new_SII'; new_SII = vertcat(new_SII{:});
else
    new_SII = {};
end
if exist('new_1')
    new_1 = new_1'; new_1 = vertcat(new_1{:});
else
    new_1 = {};
end
if exist('new_2')
    new_2 = new_2'; new_2 = vertcat(new_2{:});
else
    new_2 = {};
end
if exist('new_3')
    new_3 = new_3'; new_3 = vertcat(new_3{:});
else
    new_3 = {};
end
if exist('new_4')
    new_4 = new_4'; new_4 = vertcat(new_4{:});
else
    new_4 = {};
end
if exist('new_5')
    new_5 = new_5'; new_5 = vertcat(new_5{:});
else
    new_5 = {};
end
if exist('new_6')
    new_6 = new_6'; new_6 = vertcat(new_6{:});
else
    new_6 = {};
end
%new_HCCrejs = new_HCCrejs'; new_HCCrejs = vertcat(new_HCCrejs{:});
new_SISRDINC = new_SISRDINC'; new_SISRDINC = vertcat(new_SISRDINC{:});
new_SIDRDINC = new_SIDRDINC'; new_SIDRDINC = vertcat(new_SIDRDINC{:});




%format for the raw trace export is er-ch-time-trial
if savefile
%%SI
if any(strcmp(contr2save,'SI'))
    disp ('SI')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            disp('1')
            SI(:,:,:,:,1) = squeeze(eegpower_BN(new_si(:,1),:,:,:));
            SI(:,:,:,:,2) = squeeze(eegpower_BN(new_si(:,2),:,:,:));
        else
            disp('2')
            SI(:,:,:,:,1) = eegpower_BN(new_si(:,1),:,:,:);
            SI(:,:,:,:,2) = eegpower_BN(new_si(:,2),:,:,:);
        end
        all2all = permute (SI, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_SI_' chanNames{:} '_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end
end

%     %%DI
%     if any(strcmp(contr2save,'DI'))
%         disp ('DI')
%         if strcmp(format2save, 'pow')
%             DI(:,:,:,:,1) = squeeze(eegpower_BN(new_di(:,1), :,:,:));
%             DI(:,:,:,:,2) = squeeze(eegpower_BN(new_di(:,2), :,:,:));
%             all2all = permute (DI, [1, 5, 2, 3, 4]); 
%             filename =  [subj  '_DI_' chanNames{:} '_all2all.mat']; 
%         elseif strcmp(format2save, 'rawT')
%             DI(:,:,:,1) = oneListRooms_M (:,:,new_di(:,1));
%             DI(:,:,:,2) = oneListRooms_M (:,:,new_di(:,2));
%             all2all = permute (DI, [4, 1, 2, 3]); 
%             filename =  [subj  '_DI_' chanNames{:} '_rawT_all2all.mat']; 
%         end
%         varinfo=whos('all2all');saveopt='';
%         if (varinfo.bytes >= 2^31) saveopt ='-v7.3';end
%         save (filename, 'all2all', 'chanNames', saveopt);
%     end

%%DIN
if any(strcmp(contr2save,'DIN'))
    disp ('DIN')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            DIN(:,:,:,:,1) = squeeze(eegpower_BN(new_din(:,1), :,:,:));
            DIN(:,:,:,:,2) = squeeze(eegpower_BN(new_din(:,2), :,:,:));
        else
            DIN(:,:,:,:,1) = eegpower_BN(new_din(:,1), :,:,:);
            DIN(:,:,:,:,2) = eegpower_BN(new_din(:,2), :,:,:);
        end
        all2all = permute (DIN, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_DIN_' chanNames{:} '_all2all.mat']; 
    end
    export_batches(all2all, chanNames, '_DIN_', 200, subj); % check if needed with 0.01 resolution
end

%%DINC
if any(strcmp(contr2save,'DINC'))
    disp ('DINC')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            DINC(:,:,:,:,1) = squeeze(eegpower_BN(new_dinc(:,1), :,:,:));
            DINC(:,:,:,:,2) = squeeze(eegpower_BN(new_dinc(:,2), :,:,:));
        else
            DINC(:,:,:,:,1) = eegpower_BN(new_dinc(:,1), :,:,:);
            DINC(:,:,:,:,2) = eegpower_BN(new_dinc(:,2), :,:,:);
        end
        all2all = permute (DINC, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_DINC_' chanNames{:} '_all2all.mat']; 
    end
    export_batches(all2all, chanNames, '_DINC_', 200, subj);
end

%%DINHCC
if any(strcmp(contr2save,'DINHCC'))
    disp ('DINHCC')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            DINHCC(:,:,:,:,1) = squeeze(eegpower_BN(new_dinhcc(:,1), :,:,:));
            DINHCC(:,:,:,:,2) = squeeze(eegpower_BN(new_dinhcc(:,2), :,:,:));
        else
            DINHCC(:,:,:,:,1) = eegpower_BN(new_dinhcc(:,1), :,:,:);
            DINHCC(:,:,:,:,2) = eegpower_BN(new_dinhcc(:,2), :,:,:);
        end
        all2all = permute (DINHCC, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_DINHCC_' chanNames{:} '_all2all.mat']; 
    end
    export_batches(all2all, chanNames, '_DINHCC_', 200, subj);
end

%%SISRDINC
if any(strcmp(contr2save,'SISRDINC'))
    disp ('SISRDINC')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            SISRDINC(:,:,:,:,1) = squeeze(eegpower_BN(new_SISRDINC(:,1), :,:,:));
            SISRDINC(:,:,:,:,2) = squeeze(eegpower_BN(new_SISRDINC(:,2), :,:,:));
        else
            SISRDINC(:,:,:,:,1) = eegpower_BN(new_SISRDINC(:,1), :,:,:);
            SISRDINC(:,:,:,:,2) = eegpower_BN(new_SISRDINC(:,2), :,:,:);
        end
        all2all = permute (SISRDINC, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_SISRDINC_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')

    end
    %export_batches(all2all, chanNames, '_SISRDINC_', 200, subj);
    save (filename, 'all2all', 'chanNames');
end

%%SIDRDINC
if any(strcmp(contr2save,'SIDRDINC'))
    disp ('SIDRDINC')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            SIDRDINC(:,:,:,:,1) = squeeze(eegpower_BN(new_SIDRDINC(:,1), :,:,:));
            SIDRDINC(:,:,:,:,2) = squeeze(eegpower_BN(new_SIDRDINC(:,2), :,:,:));
        else
            SIDRDINC(:,:,:,:,1) = eegpower_BN(new_SIDRDINC(:,1), :,:,:);
            SIDRDINC(:,:,:,:,2) = eegpower_BN(new_SIDRDINC(:,2), :,:,:);
        end
        all2all = permute (SIDRDINC, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_SIDRDINC_' chanNames{:} '_all2all.mat']; 

    end
    export_batches(all2all, chanNames, '_SIDRDINC_', 1000, subj);
    %save (filename, 'all2all', 'chanNames');
end

%%SIC
if any(strcmp(contr2save,'SIC'))
    disp ('SIC')
    if strcmp(format2save, 'pow')
        SIC(:,:,:,:,1) = squeeze(eegpower_BN(new_SIC(:,1), :,:,:));
        SIC(:,:,:,:,2) = squeeze(eegpower_BN(new_SIC(:,2), :,:,:));
        all2all = permute (SIC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SIC_' chanNames{:} '_all2all.mat']; 

    end
    save (filename, 'all2all', 'chanNames');
end

%%SII
if any(strcmp(contr2save,'SII')) & ~isempty(new_SII)
    disp ('SII')
    if strcmp(format2save, 'pow')
        SII(:,:,:,:,1) = squeeze(eegpower_BN(new_SII(:,1), :,:,:));
        SII(:,:,:,:,2) = squeeze(eegpower_BN(new_SII(:,2), :,:,:));
        all2all = permute (SII, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SII_' chanNames{:} '_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');

end

%%DIC
if any(strcmp(contr2save,'DIC'))
    disp ('DIC')
    if strcmp(format2save, 'pow')
        DIC(:,:,:,:,1) = squeeze(eegpower_BN(new_DIC(:,1), :,:,:));
        DIC(:,:,:,:,2) = squeeze(eegpower_BN(new_DIC(:,2), :,:,:));
        all2all = permute (DIC, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_DIC_' chanNames{:} '_all2all.mat']; 
        save (filename, 'all2all', 'chanNames');
    end

end


%%SISP
if any(strcmp(contr2save,'SISP'))
    disp ('SISP')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            SISP(:,:,:,:,1) = squeeze(eegpower_BN(new_sisp(:,1),:,:,:));
            SISP(:,:,:,:,2) = squeeze(eegpower_BN(new_sisp(:,2),:,:,:));
        else
            SISP(:,:,:,:,1) = eegpower_BN(new_sisp(:,1),:,:,:);
            SISP(:,:,:,:,2) = eegpower_BN(new_sisp(:,2),:,:,:);
        end
        all2all = permute (SISP, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SISP_' chanNames{:} '_all2all.mat'];

    end
    save (filename, 'all2all', 'chanNames');
end

%%SIDP
if any(strcmp(contr2save,'SIDP')) & ~isempty(new_sidp)
    disp ('SIDP')
    if strcmp(format2save, 'pow')
        SIDP(:,:,:,:,1) = squeeze(eegpower_BN(new_sidp(:,1),:,:,:));
        SIDP(:,:,:,:,2) = squeeze(eegpower_BN(new_sidp(:,2),:,:,:));
        all2all = permute (SIDP, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SIDP_' chanNames{:} '_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%SIDPC
if any(strcmp(contr2save,'SIDPC')) & ~isempty(new_sidpc)
    disp ('SIDPC')
    if strcmp(format2save, 'pow')
        SIDPC(:,:,:,:,1) = squeeze(eegpower_BN(new_sidpc(:,1),:,:,:));
        SIDPC(:,:,:,:,2) = squeeze(eegpower_BN(new_sidpc(:,2),:,:,:));
        all2all = permute (SIDPC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SIDPC_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SIDPC(:,:,:,1) = oneListRooms_M (:,:,new_sidpc(:,1));
        SIDPC(:,:,:,2) = oneListRooms_M (:,:,new_sidpc(:,2));
        all2all = permute (SIDPC, [4, 1, 2, 3]); 
        filename =  [subj  '_SIDPC_' chanNames{:} '_rawT_all2all.mat'];
    end
    save (filename, 'all2all', 'chanNames');
end

%%SISR
if any(strcmp(contr2save,'SISR'))
    disp ('SISR')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            SISR(:,:,:,:,1) = squeeze(eegpower_BN(new_sisr(:,1),:,:,:));
            SISR(:,:,:,:,2) = squeeze(eegpower_BN(new_sisr(:,2),:,:,:));
        else
            SISR(:,:,:,:,1) = eegpower_BN(new_sisr(:,1),:,:,:);
            SISR(:,:,:,:,2) = eegpower_BN(new_sisr(:,2),:,:,:);
        end
        all2all = permute (SISR, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SISR_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SISR(:,:,:,1) = oneListRooms_M (:,:,new_sisr(:,1));
        SISR(:,:,:,2) = oneListRooms_M (:,:,new_sisr(:,2));
        all2all = permute (SISR, [4, 1, 2, 3]); 
        filename =  [subj  '_SISR_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%DISR
if any(strcmp(contr2save,'DISR'))
    disp ('DISR')
    if strcmp(format2save, 'pow')
        DISR(:,:,:,:,1) = squeeze(eegpower_BN(new_disr(:,1),:,:,:));
        DISR(:,:,:,:,2) = squeeze(eegpower_BN(new_disr(:,2),:,:,:));
        all2all = permute (DISR, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_DISR_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        DISR(:,:,:,1) = oneListRooms_M (:,:,new_disr(:,1));
        DISR(:,:,:,2) = oneListRooms_M (:,:,new_disr(:,2));
        all2all = permute (DISR, [4, 1, 2, 3]); 
        filename =  [subj  '_DISR_' chanNames{:} '_rawT_all2all.mat'];
    end
    save (filename, 'all2all', 'chanNames');
end

%%SIDR (incongruent items)
if any(strcmp(contr2save,'SIDR'))& ~isempty(new_sidr)
    disp ('SIDR')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            SIDR(:,:,:,:,1) = squeeze(eegpower_BN(new_sidr(:,1),:,:,:));
            SIDR(:,:,:,:,2) = squeeze(eegpower_BN(new_sidr(:,2),:,:,:));
        else
            SIDR(:,:,:,:,1) = eegpower_BN(new_sidr(:,1),:,:,:);
            SIDR(:,:,:,:,2) = eegpower_BN(new_sidr(:,2),:,:,:); 
        end

        all2all = permute (SIDR, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SIDR_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SIDR(:,:,:,1) = oneListRooms_M (:,:,new_sidr(:,1));
        SIDR(:,:,:,2) = oneListRooms_M (:,:,new_sidr(:,2));
        all2all = permute (SIDR, [4, 1, 2, 3]); 
        filename =  [subj  '_SIDR_' chanNames{:} '_rawT_all2all.mat']; 

    end
    save (filename, 'all2all', 'chanNames');
end

%%DIDR (was taking ~12 min with for loop) now ~40s
if any(strcmp(contr2save,'DIDR'))
    disp ('DIDR')
    if strcmp(format2save, 'pow')
        DIDR(:,:,:,:,1) = squeeze(eegpower_BN(new_didr(:,1),:,:,:));
        DIDR(:,:,:,:,2) = squeeze(eegpower_BN(new_didr(:,2),:,:,:));
        all2all = permute (DIDR, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_DIDR_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        DIDR(:,:,:,1) = oneListRooms_M (:,:,new_didr(:,1));
        DIDR(:,:,:,2) = oneListRooms_M (:,:,new_didr(:,2));
        all2all = permute (DIDR, [4, 1, 2, 3]); 
        filename =  [subj  '_DIDR_' chanNames{:} '_rawT_all2all.mat'];
    end
    save (filename, 'all2all', 'chanNames');
end

%%HC
if any(strcmp(contr2save,'HC')) & ~isempty(new_hc)
    disp ('HC')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 %if more than 2 channels
            HC(:,:,:,:,1) = squeeze(eegpower_BN(new_hc(:,1),:,:,:));
            HC(:,:,:,:,2) = squeeze(eegpower_BN(new_hc(:,2),:,:,:));
        else
            HC(:,:,:,:,1) = eegpower_BN(new_hc(:,1),:,:,:);
            HC(:,:,:,:,2) = eegpower_BN(new_hc(:,2),:,:,:);
        end
        all2all = permute (HC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_HC_' chanNames{:} '_all2all.mat'];
    elseif strcmp(format2save, 'rawT')
        HC(:,:,:,1) = oneListRooms_M (:,:,new_hc(:,1));
        HC(:,:,:,2) = oneListRooms_M (:,:,new_hc(:,2));
        all2all = permute (HC, [4, 1, 2, 3]); 
        filename =  [subj  '_HC_' chanNames{:} '_rawT_all2all.mat']; 

    end
    save (filename, 'all2all', 'chanNames');
end

%%LC
if any(strcmp(contr2save,'LC')) & ~isempty(new_lc)
    disp ('LC')
    if strcmp(format2save, 'pow')
        if size(eegpower_BN, 2) > 1 & size(new_lc, 1) > 1 %if more than 1 channel and more than 1 trial
            LC(:,:,:,:,1) = squeeze(eegpower_BN(new_lc(:,1),:,:,:));
            LC(:,:,:,:,2) = squeeze(eegpower_BN(new_lc(:,2),:,:,:));
        else
            LC(:,:,:,:,1) = eegpower_BN(new_lc(:,1),:,:,:);
            LC(:,:,:,:,2) = eegpower_BN(new_lc(:,2),:,:,:);
        end
        all2all = permute (LC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_LC_' chanNames{:} '_all2all.mat']; 
        save (filename, 'all2all', 'chanNames');
    elseif strcmp(format2save, 'rawT')
        LC(:,:,:,1) = oneListRooms_M (:,:,new_lc(:,1));
        LC(:,:,:,2) = oneListRooms_M (:,:,new_lc(:,2));
        all2all = permute (LC, [4, 1, 2, 3]); 
        filename =  [subj  '_LC_' chanNames{:} '_rawT_all2all.mat']; 

    end
    save (filename, 'all2all', 'chanNames');
end

if any(strcmp(contr2save,'HCC'))
    disp ('HCC')
    if strcmp(format2save, 'pow')
        HCC(:,:,:,:,1) = squeeze(eegpower_BN(new_hcc(:,1),:,:,:));
        HCC(:,:,:,:,2) = squeeze(eegpower_BN(new_hcc(:,2),:,:,:));
        all2all = permute (HCC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_HCC_' chanNames{:} '_all2all.mat'];            
        save (filename, 'all2all', 'chanNames');
    end
end

%%LCC
if any(strcmp(contr2save,'LCC')) & ~isempty(new_lcc)
    disp ('LCC')
    if strcmp(format2save, 'pow')
        LCC(:,:,:,:,1) = squeeze(eegpower_BN(new_lcc(:,1),:,:,:));
        LCC(:,:,:,:,2) = squeeze(eegpower_BN(new_lcc(:,2),:,:,:));
        all2all = permute (LCC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_LCC_' chanNames{:} '_all2all.mat']; 
        save (filename, 'all2all', 'chanNames');
    end
end


if any(strcmp(contr2save,'HCI')) & ~isempty(new_HCI)
    disp ('HCI')
    if strcmp(format2save, 'pow')
        HCI(:,:,:,:,1) = squeeze(eegpower_BN(new_HCI(:,1),:,:,:));
        HCI(:,:,:,:,2) = squeeze(eegpower_BN(new_HCI(:,2),:,:,:));
        all2all = permute (HCI, [1, 5, 2, 3, 4]); 
        filename =  [subj  '_HCI_' chanNames{:} '_all2all.mat']; 
        save (filename, 'all2all', 'chanNames');
    end
end

%%SISPC
if any(strcmp(contr2save,'SISPC'))
    disp ('SISPC')
    if strcmp(format2save, 'pow')
        SISPC(:,:,:,:,1) = squeeze(eegpower_BN(new_SISPC(:,1),:,:,:));
        SISPC(:,:,:,:,2) = squeeze(eegpower_BN(new_SISPC(:,2),:,:,:));
        all2all = permute (SISPC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SISPC_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SISPC(:,:,:,1) = oneListRooms_M (:,:,new_SISPC(:,1));
        SISPC(:,:,:,2) = oneListRooms_M (:,:,new_SISPC(:,2));
        all2all = permute (SISPC, [4, 1, 2, 3]); 
        filename =  [subj  '_SISPC_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%SIDRC
if any(strcmp(contr2save,'SIDRC')) & ~isempty(new_SIDRC)
    disp ('SIDRC')
    if strcmp(format2save, 'pow')
        SIDRC(:,:,:,:,1) = squeeze(eegpower_BN(new_SIDRC(:,1),:,:,:));
        SIDRC(:,:,:,:,2) = squeeze(eegpower_BN(new_SIDRC(:,2),:,:,:));
        all2all = permute (SIDRC, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SIDRC_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SIDRC(:,:,:,1) = oneListRooms_M (:,:,new_SIDRC(:,1));
        SIDRC(:,:,:,2) = oneListRooms_M (:,:,new_SIDRC(:,2));
        all2all = permute (SIDRC, [4, 1, 2, 3]); 
        filename =  [subj  '_SIDRC_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%1
if any(strcmp(contr2save,'SI1')) & ~isempty(new_1)
    disp ('1')
    if strcmp(format2save, 'pow')
        if size(new_1, 1) > 1
            SI1(:,:,:,:,1) = squeeze(eegpower_BN(new_1(:,1),:,:,:));
            SI1(:,:,:,:,2) = squeeze(eegpower_BN(new_1(:,2),:,:,:));
        else 
            SI1(:,:,:,:,1) = eegpower_BN(new_1(:,1),:,:,:);
            SI1(:,:,:,:,2) = eegpower_BN(new_1(:,2),:,:,:);
        end
        all2all = permute (SI1, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SI1_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SI1(:,:,:,1) = oneListRooms_M (:,:,new_1(:,1));
        SI1(:,:,:,2) = oneListRooms_M (:,:,new_1(:,2));
        all2all = permute (SI1, [4, 1, 2, 3]); 
        filename =  [subj  '_SI1_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%2
if any(strcmp(contr2save,'SI2')) & ~isempty(new_2)
    disp ('2')
    if strcmp(format2save, 'pow')
        if size(new_2, 1) > 1
            SI2(:,:,:,:,1) = squeeze(eegpower_BN(new_2(:,1),:,:,:));
            SI2(:,:,:,:,2) = squeeze(eegpower_BN(new_2(:,2),:,:,:));
        else 
            SI2(:,:,:,:,1) = eegpower_BN(new_2(:,1),:,:,:);
            SI2(:,:,:,:,2) = eegpower_BN(new_2(:,2),:,:,:);
        end
        all2all = permute (SI2, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SI2_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SI2(:,:,:,1) = oneListRooms_M (:,:,new_2(:,1));
        SI2(:,:,:,2) = oneListRooms_M (:,:,new_2(:,2));
        all2all = permute (SI2, [4, 1, 2, 3]); 
        filename =  [subj  '_SI2_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%3
if any(strcmp(contr2save,'SI3')) & ~isempty(new_3)
    disp ('3')
    if strcmp(format2save, 'pow')
        if size(new_3, 1) > 1
            SI3(:,:,:,:,1) = squeeze(eegpower_BN(new_3(:,1),:,:,:));
            SI3(:,:,:,:,2) = squeeze(eegpower_BN(new_3(:,2),:,:,:));
        else 
            SI3(:,:,:,:,1) = eegpower_BN(new_3(:,1),:,:,:);
            SI3(:,:,:,:,2) = eegpower_BN(new_3(:,2),:,:,:);
        end
        all2all = permute (SI3, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SI3_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SI3(:,:,:,1) = oneListRooms_M (:,:,new_3(:,1));
        SI3(:,:,:,2) = oneListRooms_M (:,:,new_3(:,2));
        all2all = permute (SI3, [4, 1, 2, 3]); 
        filename =  [subj  '_SI3_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%1
if any(strcmp(contr2save,'SI4')) & ~isempty(new_4)
    disp ('4')
    if strcmp(format2save, 'pow')
        if size(new_4, 1) > 1
            SI4(:,:,:,:,1) = squeeze(eegpower_BN(new_4(:,1),:,:,:));
            SI4(:,:,:,:,2) = squeeze(eegpower_BN(new_4(:,2),:,:,:));
        else 
            SI4(:,:,:,:,1) = eegpower_BN(new_4(:,1),:,:,:);
            SI4(:,:,:,:,2) = eegpower_BN(new_4(:,2),:,:,:);
        end
        all2all = permute (SI4, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SI4_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SI4(:,:,:,1) = oneListRooms_M (:,:,new_4(:,1));
        SI4(:,:,:,2) = oneListRooms_M (:,:,new_4(:,2));
        all2all = permute (SI4, [4, 1, 2, 3]); 
        filename =  [subj  '_SI4_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%1
if any(strcmp(contr2save,'SI5')) & ~isempty(new_5)
    disp ('5')
    if strcmp(format2save, 'pow')
        if size(new_5, 1) > 1
            SI5(:,:,:,:,1) = squeeze(eegpower_BN(new_5(:,1),:,:,:));
            SI5(:,:,:,:,2) = squeeze(eegpower_BN(new_5(:,2),:,:,:));
        else 
            SI5(:,:,:,:,1) = eegpower_BN(new_5(:,1),:,:,:);
            SI5(:,:,:,:,2) = eegpower_BN(new_5(:,2),:,:,:);
        end
        all2all = permute (SI5, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SI5_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SI5(:,:,:,1) = oneListRooms_M (:,:,new_5(:,1));
        SI5(:,:,:,2) = oneListRooms_M (:,:,new_5(:,2));
        all2all = permute (SI5, [4, 1, 2, 3]); 
        filename =  [subj  '_SI5_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end

%%1
if any(strcmp(contr2save,'SI6')) & ~isempty(new_6)
    disp ('6')
    if strcmp(format2save, 'pow')
        if size(new_6, 1) > 1
            SI6(:,:,:,:,1) = squeeze(eegpower_BN(new_6(:,1),:,:,:));
            SI6(:,:,:,:,2) = squeeze(eegpower_BN(new_6(:,2),:,:,:));
        else 
            SI6(:,:,:,:,1) = eegpower_BN(new_6(:,1),:,:,:);
            SI6(:,:,:,:,2) = eegpower_BN(new_6(:,2),:,:,:);
        end
        all2all = permute (SI6, [1, 5, 2, 3, 4]); size(all2all)
        filename =  [subj  '_SI6_' chanNames{:} '_all2all.mat']; 
    elseif strcmp(format2save, 'rawT')
        SI6(:,:,:,1) = oneListRooms_M (:,:,new_6(:,1));
        SI6(:,:,:,2) = oneListRooms_M (:,:,new_6(:,2));
        all2all = permute (SI6, [4, 1, 2, 3]); 
        filename =  [subj  '_SI6_' chanNames{:} '_rawT_all2all.mat']; 
    end
    save (filename, 'all2all', 'chanNames');
end





end






