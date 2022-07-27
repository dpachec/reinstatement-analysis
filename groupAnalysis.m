%% RSA GLOBAL  CREATE FOLDERS
% run in folder with _OBO_rsa files from all subj together
% first create the folders based on the conditions's names
clearvars %-except ALLEEG
sublist = dir('*rsa.mat'); sublist = {sublist.name};
fname_tmp = 'F0';
for filei=1:length(sublist)
   str = sublist{filei};
    fname_before = strsplit(str, '_');
    fname = fname_before{2};
   if ~strcmp(fname, fname_tmp)
      mkdir(fname)
      movefile(str, fname)
   else
      movefile(str, fname)
   end
   str_tmp = sublist{filei};
   str_tmp_before = strsplit(str_tmp, '_');
   fname_tmp = str_tmp_before{2};
end



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



%% plot all trials
tic
cond = SI{7};
for ti = 1:size(cond, 1)
   figure(1);
   toP = squeeze(cond(ti,:,:));
   imagesc(toP); set (gca, 'clim', [-.5 .5]); axis square;
   saveas (1, [num2str(ti) '_old.jpg']);
end
close all;
toc


%% LOAD all conditions
clearvars
tic
disp ('>>> loading conditions');
sublist = dir('*.mat');
sublist = {sublist.name};
for filei=1:length(sublist)
    load(sublist{filei});
end
toc







%% MEAN IN CLUSTER OF INTEREST
%load clustInfoReal_H ; load clustInfoReal_LTC_HC_LC;
clustInfoReal = clustInfoReal_LTC;
%clustInfoReal = clustInfoReal_H;

rsa_cond1 = cell2mat(cellfun(@mean, HC, 'un',0)); 
rsa_cond2 = cell2mat(cellfun(@mean, LC, 'un',0)); 

takeclust       = 1;
takeoverlap     = 0;
roiX  = 25:31; roiY  = 6:11; % % Hippocampus
%roiX = 25:31; roiY = 13:23; % LTC > real
pixel = 10;%H_cong_inc = 6; % LTC_HC_LC = 12 % LTC_SIDI = 16-19 // with 9 subj: H>5 ;LTC_SI>10

%without baseline
% roiX = 20:26; roiY = 9:17; % LTC > 41-75 > use pixel= 8
% pixel = 8;

mlim = 36:75;
rsa_cond1 = rsa_cond1(:,mlim,mlim);
rsa_cond2 = rsa_cond2(:,mlim,mlim);

nSubj = size(rsa_cond1, 1);
if takeclust
    ROI_cond1 = rsa_cond1(:,clustInfoReal.PixelIdxList{pixel}); 
    ROI_cond2 = rsa_cond2(:,clustInfoReal.PixelIdxList{pixel}); 
    m_roi_cond1 = mean(ROI_cond1, 2);
    m_roi_cond2 = mean(ROI_cond2, 2);
elseif takeoverlap
    ROI_cond1 = rsa_cond1(:,clustInfoOverlapping.PixelIdxList{1}); 
    ROI_cond2 = rsa_cond2(:,clustInfoOverlapping.PixelIdxList{1}); 
    m_roi_cond1_tmp = mean(ROI_cond1, 2);
    m_roi_cond2_tmp = mean(ROI_cond2, 2);
    ROI_cond1 = rsa_cond1(:,clustInfoOverlapping.PixelIdxList{2}); 
    ROI_cond2 = rsa_cond2(:,clustInfoOverlapping.PixelIdxList{2}); 
    m_roi_cond1_tmp(:,2) = mean(ROI_cond1, 2);
    m_roi_cond2_tmp(:,2) = mean(ROI_cond2, 2);
    m_roi_cond1 = mean(m_roi_cond1_tmp,2);
    m_roi_cond2 = mean(m_roi_cond2_tmp,2);
    
else
    ROI_cond1 = rsa_cond1(:,roiX, roiY); 
    ROI_cond2 = rsa_cond2(:,roiX, roiY); 
    ROI_cond1 = reshape (ROI_cond1, [nSubj size(roiX,2) * size(roiY, 2)]);
    ROI_cond2 = reshape (ROI_cond2, [nSubj size(roiX,2) * size(roiY, 2)]);
    m_roi_cond1 = mean(ROI_cond1, 2);
    m_roi_cond2 = mean(ROI_cond2, 2);
end

test = squeeze(mean(rsa_cond1, 1));
test(1:end) = 0;
%test(roiX,roiY,:) = .5;
if takeclust
    test(clustInfoReal.PixelIdxList{pixel}) = 1;
elseif takeoverlap
    test(clustInfoOverlapping.PixelIdxList{1}) = 1;
end
test = flipud(test);
figure(1);imagesc(test);axis square;hold on;
%figure(1);contourf(test, 40, 'linecolor', 'none');axis square;hold on;
mfL = 5; bins = 40;
plot([mfL+0.5 mfL+0.5],get(gca,'ylim'),'w', 'LineWidth', 2); 
plot(get(gca,'xlim'), [bins - (mfL-0.5) bins - (mfL-0.5)],'w', 'LineWidth', 2);
plot(get(gca,'xlim'), [bins+.5 .5],'w', 'LineWidth', 2); % diagonal
set (gca, 'clim', [0 1]);
%remove labels
axesHandles = findall(0, 'type', 'axes');
for i=1:length(axesHandles)
   set (axesHandles(i), 'visible', 'off'); % remove labels 
end


[h p ci t] = ttest (m_roi_cond1);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set (gcf, 'Position', [200 200 400 450]);
export_fig(2, 'plotting1.png','-transparent', '-r300');
close all;   

%% PLOT 2 bar
clear data.data
%ylim = [-.035 .08];
ylim = [-.125 .15]; %for the hc-lc comparison

data.data = [m_roi_cond1 m_roi_cond2];
figure(2); set(gcf,'Position', [0 0 560 500]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(7);
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
hb = plot ([1 2 ], data.data(:,1:2)); hold on; % > lines
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
%set(hb,'linestyle','none', 'lineWidth', 2);
%set(hb, 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LC  '}, ...
    'FontSize', 32, 'linew',3, 'xlim', [0.4 2.6], 'ylim', ylim );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);
xticklabel_rotate;

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

export_fig(2, '_currentBarPlot_m2C.png','-transparent', '-r300');
close all;   

%% PLOT 1 bar
clear data.data
ylim = [-.01 .1];
%ylim = [-.125 .15]; %for the hc-lc comparison

data.data = [m_roi_cond1];
figure(2); set(gcf,'Position', [0 0 560 500]); 
mean_S = mean(data.data, 1);
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
hb = plot ([1 ], data.data(:,1)); hold on; % > lines
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',50);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2); 
set(gca,'XTick',[1],'XTickLabel',{'HC  '}, ...
    'FontSize', 32, 'linew',3, 'xlim', [0 2], 'ylim', ylim );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

export_fig(2, 'oneBar.png','-transparent', '-r300');
close all;   




%% PLOT 6 bar
ylim = [-.035 .075];
%ylim = [-.1 .13]; %for the hc-lc comparison
load all_conditions_in_H_SISP-SIDR_cluster

data.data = [SISP SIDR SI DI SR DR];
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(7);
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
hb = plot ([1 2 3 4 5 6 ], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
%set(hb, 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'Congruent  ', 'Incongruent  ' ...
    'Same Item  ', 'Different Item  ', 'Same Room  ', 'Different Room  '}, ...
    'FontSize', 25, 'linew',2, 'xlim', [0 7], 'ylim', ylim );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);
xticklabel_rotate;

[h p ci t] = ttest (data.data(:,6));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);


set (gcf, 'Position', [200 200 500 550]);
export_fig(2, '_currentBarPlot_m2C.png','-transparent', '-r300');
close all;   



%% LOAD DATA and exclude subjects from the saved contrasts
clearvars; 
load SI_H; load SI_LTC; load SISP_H; load SISP_LTC; load SIDR_H; 
load SIDR_LTC; load clustInfoReal_H; load clustInfoReal_LTC_HC_LC; 
s = [2];
HC{s} = []; HC = HC(~cellfun('isempty',HC)); 
SI_H{s} = []; SI_H = SI_H(~cellfun('isempty',SI_H)); 
SI_LTC{s} = []; SI_LTC = SI_LTC(~cellfun('isempty',SI_LTC)); 
SISP_H{s} = []; SISP_H = SISP_H(~cellfun('isempty',SISP_H)); 
SISP_LTC{s} = []; SISP_LTC = SISP_LTC(~cellfun('isempty',SISP_LTC)); 
SIDR_H{s} = []; SIDR_H = SIDR_H(~cellfun('isempty',SIDR_H)); 
SIDR_LTC{s} = []; SIDR_LTC = SIDR_LTC(~cellfun('isempty',SIDR_LTC)); 

%% COORDINATED REINSTATEMENT IN SPECIFIC CLUSTERS
%clearvars;
load clustInfoReal_H ; load clustInfoReal_LTC_HC_LC;
load clustInfoReal_LTC_SI_DI;

saveimg = 1;
takeClusters = 0;
clear all_trials_RSA_H all_trials_RSA_LTC;
all_trials_RSA_H = SIDR_H;
all_trials_RSA_LTC = SIDR_LTC;

%40
pixel_LTC = 12;% HC-LC_in_LTC = 12;pixel = 17; LTC_SIDI = 17-20
mlim = 36:75;

%%time matched 
pixel_H = 6; %H_cong_inc = 6; %1 for matched at encoding %
roiX_H  = 25:31; roiY_H  = 6:11; % 
roiX_LTC = 25:31; roiY_LTC  = 13:23; % 

%%most significant 
% pixel_H = 14; %use biggestSmallest H cluster. H_cong_inc = 6; %1 for matched at encoding %
% roiX_H  = 23:28; roiY_H  = 7:11; % 
% roiX_LTC = 26:31; roiY_LTC  = 14:22; % 

clear perTrialRSA_H perTrialRSA_LTC subjt_H subjt_LTC;
for subji = 1:length(all_trials_RSA_H)
    allTRSA_H = all_trials_RSA_H{subji}; 
    allTRSA_H = allTRSA_H (:, mlim, mlim);
    allTRSA_LTC = all_trials_RSA_LTC{subji};
    allTRSA_LTC = allTRSA_LTC (:, mlim, mlim);
    for ti = 1:size(allTRSA_H, 1)
        tMat = allTRSA_H(ti, roiX_H, roiY_H); %square
        if (takeClusters)
            tMat = allTRSA_H(ti, clustInfoReal_H.PixelIdxList{pixel_H}); %cluster
        end
        t = tMat(:);   t = mean (t); 
        subjt_H(ti,:) = t;
        perTrialRSA_H{subji} = subjt_H'; 
    end
    for ti = 1:size(allTRSA_LTC, 1)
        t = allTRSA_LTC(ti, roiX_LTC, roiY_LTC);  %square
        if (takeClusters)
            t = allTRSA_LTC(ti, clustInfoReal_LTC.PixelIdxList{pixel_LTC}); %cluster
        end
        t = t(:);   t = mean (t); %cluster
        subjt_LTC(ti,:) = t;
        perTrialRSA_LTC{subji} = subjt_LTC'; 
    end
    
    clear subjt_H subjt_LTC;
end


perTrialRSA_H = perTrialRSA_H'; perTrialRSA_LTC = perTrialRSA_LTC';

%%CORRELATE RSA per trial in H and LTC
clear rsa;
for subji = 1:length(perTrialRSA_H)
    x = perTrialRSA_H{subji}(:);
    y = perTrialRSA_LTC{subji}(:);
    [r,p] = corr(x, y, 'Type','s'); %Kendall
    rsa(subji,1) = r; rsa(subji,2) = p;
    meanX(subji) = mean(x);
    meanY(subji) = mean(y);
end

rsaZ = atanh(squeeze(rsa));
rsaZ(:,1)'
% t-test
[h p ci t] = ttest (rsaZ(:,1)); %1 means Rhos 
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

clear test;
test = squeeze(mean(allTRSA_LTC, 1));
test(1:end) = 0.2;
%squares
test(roiX_H,roiY_H,:) = .55;
test(roiX_LTC,roiY_LTC,:) = .55;

%pixels
test(clustInfoReal_H.PixelIdxList{pixel_H}) = 1; 
%test(clustInfoReal_LTC.PixelIdxList{pixel_LTC}) = 1;
test = flipud(test); 
figure(1);imagesc(test);axis square;hold on;
%title(['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
%remove labels
axesHandles = findall(0, 'type', 'axes');
for i=1:length(axesHandles)
   set (axesHandles(i), 'visible', 'off'); % remove labels 
end



if saveimg
titFig = 'current_fig';
export_fig (titFig, '-transparent', '-r200');
close all;
end

%% plot one bar
ylim = [-.1 .25];
xlim = [0 2];

data.data = [rsaZ(:,1)];
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',45);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'All items'}, ...
    'FontSize', 30, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);
%xticklabel_rotate;

set (gcf, 'Position', [200 200 450 450]);
export_fig(2, '_currentBarPlot.png','-transparent', '-r300');
close all; 


%% check for differences in coordination for cong and inc trials
load coord_SIDR; load coord_SISP;
[h p ci t] = ttest(SIDR(:, 1), SISP(:,1))


%% SLIDING 40 (hippocampal seed)
tic
clearvars -except SI_H SI_LTC SISP_H SISP_LTC clustInfoReal_H clustInfoReal_LTC roiX_H roiY_H;
%roiX_H = 24:29; roiY_H  = 6:11; % ROIX has to be defined in the previous
%cell

takeClust = 0;
pixel2use = 6;
mlim = 36:75;
all_trials_RSA_H = SISP_H;
all_trials_RSA_LTC = SISP_LTC;
nSubj = length(all_trials_RSA_H);
tic

xlen = 33;
ylen = 31;

for subji = 1:nSubj
        allTRSA_H = all_trials_RSA_H{subji}; 
        allTRSA_H = allTRSA_H (:, mlim, mlim);
        allTRSA_LTC = all_trials_RSA_LTC{subji};
        allTRSA_LTC = allTRSA_LTC (:, mlim, mlim);
        nTrials = size(all_trials_RSA_H{subji},1);
        perTrialRSA_H = zeros (1, nTrials);
        perTrialRSA_LTC = zeros (nTrials, xlen,ylen);
        
        if takeClust
            tH = allTRSA_H(:, clustInfoReal_H.PixelIdxList{pixel2use}); %square
            tvecH = reshape(tH, [nTrials, length(clustInfoReal_H.PixelIdxList{pixel2use})]); %I need the mean for each trial
        else
            tH = allTRSA_H(:, roiX_H, roiY_H); %square
            tvecH = reshape(tH, [nTrials, size(roiX_H,2)*size(roiY_H,2)]); %I need the mean for each trial
        end
        tM_H = mean(tvecH, 2); %> tm mean in ROI
        perTrialRSA_H = tM_H; 
        
    for xi = 1:xlen
        fprintf('%d %s', xi, ' '); if (rem(xi,25) == 0) fprintf('\n'); end
        xV = 7; %change xi to allow reshape
        yV = 9; %change yi to allow reshape
        for yi = 1:ylen
            roiX_LTC = xi:xi+xV; roiY_LTC =  yi:yi+yV;
            tLTC = allTRSA_LTC(:, roiX_LTC, roiY_LTC); %square
            tvecLTC = reshape(tLTC, [nTrials, (xV +1) * (yV + 1)]); %I need the mean for each trial
            tM_LTC = mean(tvecLTC, 2); %> tm mean in ROI
            [r,p] = corr(tM_H, tM_LTC, 'Type','Spearman'); 
            
            perTrialRSA_H_LTC_R(:, xi, yi) = r;
            perTrialRSA_H_LTC_P(:, xi, yi) = p;
         
            %check it in the x axis with figure
%             if subji == 1 & xi < 3 & yi > 28%~mod(yi,2) 
%                 figure();
%                 tmp = squeeze(mean(allTRSA_H));
%                 tmp (roiX_LTC, roiY_LTC) = 100;
%                 tmp = flipud(tmp);
%                 imagesc (tmp); 
%             end
            %check it in the y axis with figure
%             if subji == 1 & xi > 31 & yi < 3%~mod(yi,2) 
%                 disp('hola');
%                 figure();
%                 tmp = squeeze(mean(allTRSA_H));
%                 tmp (roiX_LTC, roiY_LTC) = 100;
%                 tmp = flipud(tmp);
%                 imagesc (tmp); 
%             end
           
        end
    end
    fprintf('\n');
    allSubjRSA_H_LTC(subji,1, :,:) = perTrialRSA_H_LTC_R;
    allSubjRSA_H_LTC(subji,2, :,:) = perTrialRSA_H_LTC_P;
end

rsa = squeeze(allSubjRSA_H_LTC (:,1,:,:)); %% 1:rho, 2: p-values
rsaZ = atanh(squeeze(rsa));
[h p ci t] = ttest (rsaZ, 0);
hmap = squeeze(h); pmap = squeeze(p); tmap = squeeze(t.tstat);

filename = ['sliding_real.mat'];
save (filename, 'pmap', 'tmap', 'hmap'); 


%% SLIDING 40 (LTC seed)
tic
clearvars -except SI_H SI_LTC clustInfoReal_H clustInfoReal_LTC ...
    roiX_LTC roiY_LTC SISP_H SISP_LTC SIDR_H SIDR_LTC;
%ROILTC has to be defined in previous cell
%cell

takeClust = 0;
pixel2use = 6;
plotsquares2check = 0;
mlim = 36:75;
all_trials_RSA_H = SISP_H;
all_trials_RSA_LTC = SISP_LTC;
nSubj = length(all_trials_RSA_H);
tic

xlen = 34;
ylen = 35;

for subji = 1:nSubj
        allTRSA_H = all_trials_RSA_H{subji}; 
        allTRSA_H = allTRSA_H (:, mlim, mlim);
        allTRSA_LTC = all_trials_RSA_LTC{subji};
        allTRSA_LTC = allTRSA_LTC (:, mlim, mlim);
        nTrials = size(all_trials_RSA_H{subji},1);
        perTrialRSA_H = zeros (1, nTrials);
        perTrialRSA_LTC = zeros (nTrials, xlen,ylen);
        
        %define seed cluster
        if takeClust
            tLTC = allTRSA_LTC(:, clustInfoReal_LTC.PixelIdxList{pixel2use}); %square
            tvecLTC = reshape(tLTC, [nTrials, length(clustInfoReal_LTC.PixelIdxList{pixel2use})]); %I need the mean for each trial
        else
            tLTC = allTRSA_LTC(:, roiX_LTC, roiY_LTC); %square
            tvecLTC = reshape(tLTC, [nTrials, size(roiX_LTC,2)*size(roiY_LTC,2)]); %I need the mean for each trial
        end
        tM_LTC = mean(tvecLTC, 2); %> tm mean in ROI
        perTrialRSA_LTC = tM_LTC; 
        
        %roiX_H  = 25:31; roiY_H  = 6:11; % 
    for xi = 1:xlen
        fprintf('%d %s', xi, ' '); if (rem(xi,25) == 0) fprintf('\n'); end
        xV = 7; %change xi to allow reshape
        yV = 6; %change yi to allow reshape
        for yi = 1:ylen
            roiX_H = xi:xi+xV-1; roiY_H =  yi:yi+yV-1;
            tH = allTRSA_H(:, roiX_H, roiY_H); %square
            tvecH = reshape(tH, [nTrials, (xV) * (yV)]); %I need the mean for each trial
            tM_H = mean(tvecH, 2); %> tm mean in ROI
            
            %calculate correlation
            [r,p] = corr(tM_H, tM_LTC, 'Type','s'); 
            perTrialRSA_H_LTC_R(:, xi, yi) = r;
            perTrialRSA_H_LTC_P(:, xi, yi) = p;
         
            if plotsquares2check
                %%%check it in the x axis with figure
                if subji == 1 & xi < 3 & yi > 31 
                    figure();
                    tmp = squeeze(mean(allTRSA_H));
                    tmp (roiX_H, roiY_H) = 100;
                    tmp = flipud(tmp);
                    imagesc (tmp);  axis square;
                end
                %%%check it in the y axis with figure
                if subji == 1 & xi > 31 & yi < 3 
                    figure();
                    tmp = squeeze(mean(allTRSA_H));
                    tmp (roiX_H, roiY_H) = 100;
                    tmp = flipud(tmp);
                    imagesc (tmp); axis square;                    
                end
            end
        end
    end
    fprintf('\n');
    allSubjRSA_H_LTC(subji,1, :,:) = perTrialRSA_H_LTC_R;
    allSubjRSA_H_LTC(subji,2, :,:) = perTrialRSA_H_LTC_P;
end

rsa = squeeze(allSubjRSA_H_LTC (:,1,:,:)); %% 1:rho, 2: p-values
rsaZ = atanh(squeeze(rsa));
[h p ci t] = ttest (rsaZ, 0);
hmap = squeeze(h); pmap = squeeze(p); tmap = squeeze(t.tstat);

filename = ['sliding_real.mat'];
save (filename, 'pmap', 'tmap', 'hmap'); 

%% plot real 
% tmap (end:40, :) = 0; tmap (:, end:40) = 0;
% hmap (end:40, :) = 0; hmap (:, end:40) = 0;
tmap = tmap (3:end-2, 3:end-3);hmap = hmap (3:end-2, 3:end-3);
%hmap (end:40, :) = 0; hmap (:, end:40) = 0;
figure(3); imagesc (flipud(tmap)); hold on; colormap jet; axis square
contour(flipud(hmap), 1, 'LineWidth',2, 'linecolor', 'k');%colorbar;
set (gca, 'clim', [-3 4]); %
%remove labels
axesHandles = findall(0, 'type', 'axes');
    for i=1:length(axesHandles)
       set (axesHandles(i), 'visible', 'off'); % remove labels 
    end
export_fig(3, '_coord_reinst.png','-transparent', '-r300');
%close all;

%% plot grid
grid = rand(40, 40);
figure(2); imagesc (grid); hold on; axis square;
axesHandles = findall(0, 'type', 'axes');
    for i=1:length(axesHandles)
       set (axesHandles(i), 'visible', 'off'); % remove labels 
    end
export_fig(2, 'grid.png','-transparent', '-r300');
close all; 


%% Exctract real max cluster size
clustinfoReal = bwconncomp(hmap); 
numPixReal = cellfun(@numel,clustinfoReal.PixelIdxList);
[max_clust_infoR maxiR] = max([ cellfun(@numel,clustinfoReal.PixelIdxList) ]); % the zero accounts for empty maps
      
%%get all sums of t's
clear all_clust_tsum_real;
for ci = 1:length(clustinfoReal.PixelIdxList)
    all_clust_tsum_real(ci) = sum (tmap(clustinfoReal.PixelIdxList{ci}));
end
all_clust_tsum_real = all_clust_tsum_real';
all_clust_tsum_real(:,2) =  numPixReal;
    
max_clust_sum_real = sum (tmap(clustinfoReal.PixelIdxList{maxiR}));
max_clust_avg_real = sum (tmap(clustinfoReal.PixelIdxList{maxiR})) / ...
            length (tmap(clustinfoReal.PixelIdxList{maxiR}));


%% SLIDING 40 Permutation (both methods, LTC seed)
clearvars -except SI_H SI_LTC clustInfoReal_H clustInfoReal_LTC roiX_LTC roiY_LTC ... 
    max_clust_sum_real max_clust_avg_real numPixReal;
%ROI info has to be defined two cells above

tic
nPerm = 500;
permMethod = 'trials'; %'chans' or 'trials'
takeClust = 0;

mlim = 36:75;
all_trials_RSA_H = SI_H;
all_trials_RSA_LTC = SI_LTC;
nSubj = length(all_trials_RSA_H);
tic

xlen = 34;
ylen = 35;

fprintf('\n'); fprintf('Permutation:        '); %check this later
for permi = 1:nPerm
if (permi < 10) fprintf('\b\b\b'); elseif (permi < 100) fprintf('\b\b\b\b'); 
else fprintf('\b\b\b\b\b');  end
fprintf('%d %s', permi, ' '); 

for subji = 1:nSubj
    
    %define method to shuffle labels
    if strcmp (permMethod, 'chans')
        allTRSA_H = all_trials_RSA_H{subji}; 
        allTRSA_LTC = all_trials_RSA_LTC{subji};
        all_H_LTC = cat (1, allTRSA_H, allTRSA_LTC);
        all_H_LTC_fake = all_H_LTC(randperm(size(all_H_LTC,1)), :, :);
        allTRSA_H = all_H_LTC_fake(1:round(size(all_H_LTC,1)/2), :, :);
        allTRSA_LTC = all_H_LTC_fake(round(size(all_H_LTC,1)/2)+1:end, :, :);
        nTrials = size (allTRSA_H, 1);
    end
    
    if strcmp (permMethod, 'trials')
        allTRSA_H = all_trials_RSA_H{subji}; 
        allTRSA_LTC = all_trials_RSA_LTC{subji};
        allTRSA_H = allTRSA_H(randperm(size(allTRSA_H,1)), :, :);
        allTRSA_LTC = allTRSA_LTC(randperm(size(allTRSA_LTC,1)), :, :);
        nTrials = size (allTRSA_H, 1);
    end
    
        
    %define seed cluster
    if takeClust
        tLTC = allTRSA_LTC(:, clustInfoReal_LTC.PixelIdxList{pixel2use}); %square
        tvecLTC = reshape(tLTC, [nTrials, length(clustInfoReal_LTC.PixelIdxList{pixel2use})]); %I need the mean for each trial
    else
        tLTC = allTRSA_LTC(:, roiX_LTC, roiY_LTC); %square
        tvecLTC = reshape(tLTC, [nTrials, size(roiX_LTC,2)*size(roiY_LTC,2)]); %I need the mean for each trial
    end
    tM_LTC = mean(tvecLTC, 2); %> tm mean in ROI
    perTrialRSA_LTC = tM_LTC; 
        
    for xi = 1:xlen
        %fprintf('%d %s', xi, ' '); if (rem(xi,25) == 0) fprintf('\n'); end
        xV = 7; %change xi to allow reshape
        yV = 6; %change yi to allow reshape
        for yi = 1:ylen
            roiX_H = xi:xi+xV-1; roiY_H =  yi:yi+yV-1;
            tH = allTRSA_H(:, roiX_H, roiY_H); %square
            tvecH = reshape(tH, [nTrials, (xV) * (yV)]); %I need the mean for each trial
            tM_H = mean(tvecH, 2); %> tm mean in ROI
            
            %calculate correlation
            [r,p] = corr(tM_H, tM_LTC, 'Type','Spearman'); 
            perTrialRSA_H_LTC_R(:, xi, yi) = r;
            perTrialRSA_H_LTC_P(:, xi, yi) = p;
          
        end
    end
    allSubjRSA_H_LTC(subji,1, :,:) = perTrialRSA_H_LTC_R;
    allSubjRSA_H_LTC(subji,2, :,:) = perTrialRSA_H_LTC_P;
end

rsa = squeeze(allSubjRSA_H_LTC (:,1,:,:)); %% 1:rho, 2: p-values
rsaZ = atanh(squeeze(rsa));
[h p ci t] = ttest (rsaZ, 0);
hmap = squeeze(h); pmap = squeeze(p); tmap = squeeze(t.tstat);

permRsas_p (permi, :, :) = pmap; 
permRsas_t (permi, :, :) = tmap; 
permRsas_h (permi, :, :) = hmap; 

end
fprintf('\n');

toc

filename = [num2str(nPerm) 'p_coord_reinst.mat'];
save (filename, 'permRsas_p', 'permRsas_t', 'permRsas_h'); 




%% SLIDING 40 Permutation (both methods, Hippocampal seed)
clearvars -except SI_H SI_LTC clustInfoReal_H clustInfoReal_LTC roiX_H roiY_H ... 
    max_clust_sum_real max_clust_avg_real numPixReal;
%ROI info has to be defined two cells above

tic
nperm = 1000;
permMethod = 'chans'; %'chans' or 'trials'
takeClust = 0;

mlim = 36:75;
all_trials_RSA_H = SI_H;
all_trials_RSA_LTC = SI_LTC;
nSubj = length(all_trials_RSA_H);
tic

xlen = 34;
ylen = 32;

fprintf('\n'); fprintf('Permutation:        '); %check this later
for permi = 1:nperm
if (permi < 10) fprintf('\b\b\b'); elseif (permi < 100) fprintf('\b\b\b\b'); 
else fprintf('\b\b\b\b\b');  end
fprintf('%d %s', permi, ' '); 

for subji = 1:nSubj
    if strcmp (permMethod, 'chans')
        allTRSA_H = all_trials_RSA_H{subji}; 
        allTRSA_LTC = all_trials_RSA_LTC{subji};
        all_H_LTC = cat (1, allTRSA_H, allTRSA_LTC);
        all_H_LTC_fake = all_H_LTC(randperm(size(all_H_LTC,1)), :, :);
        allTRSA_H = all_H_LTC_fake(1:round(size(all_H_LTC,1)/2), :, :);
        allTRSA_LTC = all_H_LTC_fake(round(size(all_H_LTC,1)/2)+1:end, :, :);
        nTrials = size (allTRSA_H, 1);
    end
    
    if strcmp (permMethod, 'trials')
        allTRSA_H = all_trials_RSA_H{subji}; 
        allTRSA_LTC = all_trials_RSA_LTC{subji};
        allTRSA_H = allTRSA_H(randperm(size(allTRSA_H,1)), :, :);
        allTRSA_LTC = allTRSA_LTC(randperm(size(allTRSA_LTC,1)), :, :);
        nTrials = size (allTRSA_H, 1);
    end
    
        
    if takeClust
        tH = allTRSA_H(:, clustinfoReal_H.PixelIdxList{4}); %square
        tvecH = reshape(tH, [nTrials, length(clustinfoReal_H.PixelIdxList{4})]); %I need the mean for each trial
    else
        tH = allTRSA_H(:, roiX_H, roiY_H); %square
        tvecH = reshape(tH, [nTrials, size(roiX_H,2)*size(roiY_H,2)]); %I need the mean for each trial
    end
        tM_H = mean(tvecH, 2); %> tm mean in ROI
        perTrialRSA_H = tM_H; 
        
    for xi = 1:xlen
        %fprintf('%d %s', xi, ' '); if (rem(xi,25) == 0) fprintf('\n'); end
        xV = 7; %change xi to allow reshape
        yV = 9; %change yi to allow reshape
        for yi = 1:ylen
            roiX_LTC = xi:xi+xV; roiY_LTC =  yi:yi+yV;
            tLTC = allTRSA_LTC(:, roiX_LTC, roiY_LTC); %square
            tvecLTC = reshape(tLTC, [nTrials, (xV +1) * (yV + 1)]); %I need the mean for each trial
            tM_LTC = mean(tvecLTC, 2); %> tm mean in ROI
            [r,p] = corr(tM_H, tM_LTC, 'Type','Spearman'); 
            
            perTrialRSA_H_LTC_R(:, xi, yi) = r;
            perTrialRSA_H_LTC_P(:, xi, yi) = p;
          
        end
    end
    allSubjRSA_H_LTC(subji,1, :,:) = perTrialRSA_H_LTC_R;
    allSubjRSA_H_LTC(subji,2, :,:) = perTrialRSA_H_LTC_P;
end

rsa = squeeze(allSubjRSA_H_LTC (:,1,:,:)); %% 1:rho, 2: p-values
rsaZ = atanh(squeeze(rsa));
[h p ci t] = ttest (rsaZ, 0);
hmap = squeeze(h); pmap = squeeze(p); tmap = squeeze(t.tstat);

permRsas_p (permi, :, :) = pmap; 
permRsas_t (permi, :, :) = tmap; 
permRsas_h (permi, :, :) = hmap; 

end
fprintf('\n');
filename = ['coord_reinst_perm.mat'];
save (filename, 'permRsas_p', 'permRsas_t', 'permRsas_h'); 




%% plot random sample of permutations
x = randsample(nPerm, 20);
%x = 1:20;
figure();
for rp1 = 1:20
    subplot (4, 5,rp1)
    pmap_p = squeeze (permRsas_p (x(rp1), :, :));
    tmap_p = squeeze (permRsas_t (x(rp1), :, :));
    hmap_p = squeeze (permRsas_h (x(rp1), :, :));
    imagesc (flipud(tmap_p)); hold on; colormap jet;
    contour(flipud(hmap_p), 1, 'LineWidth',2, 'linecolor', 'k');axis square; 
    set (gca, 'clim', [-3 7]); hold on;
end

%% cut permutation space
permRsas_h = permRsas_h (:, 3:end-2, 3:end-3);
permRsas_p = permRsas_p (:, 3:end-2, 3:end-3);
permRsas_t = permRsas_t (:, 3:end-2, 3:end-3);



%% compare real and permuted values

%%extract cluster size for each permutation and rank
pval = 0.05;
for permi = 1:size(permRsas_h,1)
    % get number of elements in largest supra-threshold cluster
    sigMT_permi = permRsas_t(permi, :, :);
    sigMH_permi = permRsas_h(permi, :, :);
    clustinfo = bwconncomp(sigMH_permi);
    [numPixPermi(permi) maxi] = max([0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
    if numPixPermi(permi) > 0
        max_clust_sum(permi) = sum (sigMT_permi(clustinfo.PixelIdxList{maxi-1}));
        max_clust_avg(permi) = max_clust_sum(permi) / ...
            length (sigMT_permi(clustinfo.PixelIdxList{maxi-1}));
    end
end

numPixReal = numPixReal';
numPixPermi = numPixPermi';
max_clust_sum = max_clust_sum'; 
max_clust_avg = max_clust_avg'; 

%%threshold to the sum of the significant cluster of t values
lower_threshold = prctile(max_clust_sum,    (pval*100)/2);
%max_clust_sum_perm = prctile(max_clust_sum,100-((pval*100)/2)); % /2 ?
max_clust_sum_perm = prctile(max_clust_sum,100-((pval*100)));
max_clust_avg_perm = prctile(max_clust_avg,100-((pval*100)));

numPixThres = prctile(abs(numPixPermi),100-pval*100);
whichclusters2remove = find(numPixReal<numPixThres);
sigMH_thres = hmap;

%remove small clusters
% for i=1:length(whichclusters2remove)
%     sigMH_thres(clustinfoReal.PixelIdxList{whichclusters2remove(i)})=0;
% end



%% get rank (for stats)
max_clust_sum_ranked = max_clust_sum;
max_clust_sum_ranked(:,2) =  tiedrank (max_clust_sum);

[c index] = min(abs(max_clust_sum_ranked(:,1)-max_clust_sum_real));
%[c index] = min(abs(max_clust_sum_ranked(:,1)- 1.163954240917805e+02 ));
1000 - max_clust_sum_ranked(index, 2)
%x = find(max_clust_sum_ranked(:,2) == index)


%% import Phase sync
%can be used to import global and local RSA
%clear, close all
clearvars 
tic

sublist = dir('*AVG.mat');%use here PLI
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    load(sublist{subji});
    %chan2plot(subji)
    size(ispc)
    blne = [36:36]; %this blne works= [32:38];
    %ispcB = 100* bsxfun(@rdivide, ispc - mean(ispc(:,:,blne(1):blne(end)), 3), mean(ispc(:,:,blne(1):blne(end)), 3)); %pick first trial just to check
    ispcB = ispc - mean(ispc(:,:,blne(1):blne(end)), 3);
    all{subji,1} = ispc; %when more than 1 electrode
    %all{subji,1} = rsaZ;
    
end

if exist ('timeBins')
    timeBins1 = timeBins (:, [1, end]);
end


cd .. % goes up one directory
filename = [sublist{subji}(4:end-8) '_pha'];
%filename = [sublist{subji}(5:end-12)];
eval([filename '= all;']);
save (filename, filename, '-v7.3');
%save (filename, filename);

toc

%% remove sub2
clearvars;
load SI_pha; load SI_rsa_LTC
SI_pha{2} = [];
SI_rsa{2} = [];
SI_pha = SI_pha(~cellfun('isempty',SI_pha));
SI_rsa = SI_rsa(~cellfun('isempty',SI_rsa));
disp ('subj removed');
%% plot mean
all = SI_pha;
m_all = cell2mat(cellfun(@mean, all, 'un',0)); 
[h p ci t] = ttest (m_all, 0, 'alpha', 0.05); 
h = squeeze(h);
t = squeeze(t.tstat);
%44F
%min_freq1 = 1;max_freq1 = 29;num_frex1 = 29;min_freq2 = 30;max_freq2 = 100;num_frex2 = 15;baseline_window = [ -500 0]; frex1 = linspace(min_freq1,max_freq1,num_frex1);frex2 = linspace(min_freq2,max_freq2,num_frex2);frex = [frex1 frex2]; num_frex = size(frex, 2);
%100F
frex = 1:100;

grand_mean_all = squeeze(mean(m_all));
%times = -3.75:.1:3.75;
times = 1:76;

figure();
subplot (211);
contourf(times, frex, grand_mean_all, 40, 'linecolor', 'none'); axis square;hold on
%plot([0 0],get(gca,'ylim'),'k');hold on;
plot([41 41],get(gca,'ylim'),'k');hold on;
%set(gca, 'xlim', [-.5 3], 'clim', [-.025 .025]); colorbar; %colorbar, 'ylim', [1 100],
set(gca, 'xlim', [35 75], 'clim', [-.025 .025]); colorbar %, 'ylim', [1 100],

subplot(212)
contourf(times, frex, t, 40, 'linecolor', 'none'); axis square;hold on
contour (times, frex, h, 1, 'lineWidth', 1, 'linecolor', 'k');
%plot([0 0],get(gca,'ylim'),'k');hold on;
plot([41 41],get(gca,'ylim'),'k');hold on;
%set(gca, 'xlim', [-.5 3], 'clim', [-3 3]); colorbar %, 'ylim', [1 100],
set(gca, 'xlim', [35 75], 'clim', [-3 3]); colorbar %, 'ylim', [1 100],

export_fig(1, '_phase_sync.png','-transparent', '-r300');

%% get mean in cluster for each subject
for subji = 1:7
   meanInClust(subji, :) =  squeeze (mean (mean (SI_pha{subji}(:,70:83,41), 2)));
    
end
[h p] = ttest(meanInClust)
figure();
boxplot(meanInClust);
xlabel ('All trials'); ylabel ('Average phase diff');
set (gca, 'FontSize', 18);



%%
t = grand_mean_all(70:83, 41)
[h p] = ttest(t);
boxplot(t)

%% calculate per trial RSA
clearvars;
load SI_rsa_LTC;
cond2use_rsa    = SI_rsa;

nSubj = length(cond2use_rsa);
clear perTrialRSA;
mlim = 36:75
roiX = 25:32; roiY = 13:23; % LTC
for subji = 1:nSubj
    rsa_map = cond2use_rsa{subji}(:,mlim, mlim);
    nTrials = size(rsa_map, 1);
    clear rsa2cmp_tmp;
    for triali = 1:nTrials
        rsaTmp = rsa_map(triali, roiX, roiY); rsaTmp = rsaTmp(:);
        rsa2cmp_tmp(triali) = mean(rsaTmp);        
    end
    perTrialRSA{subji} = rsa2cmp_tmp;
end
perTrialRSA = perTrialRSA';


%% correlation with RSA
clearvars;
tic;
load perTrialRSA_LTC_SI; load SI_pha;

takeClust       = 0;
cond2use_pha    = SI_pha;
runPerm         = 0;
nPerm           = 1000;

%exclude s2
cond2use_pha{2} = []; perTrialRSA{2} = [];
cond2use_pha = cond2use_pha(~cellfun('isempty',cond2use_pha));
perTrialRSA = perTrialRSA(~cellfun('isempty',perTrialRSA));
disp ('subj 2 removed');

%define mlim for phase map
mlim = 41:58; %use 75 for all time (3.5s) or 58 until the end of the LTC ROI
%mean rsa in cluster
nSubj = length(cond2use_pha);

% apply limits to phase data
for subji = 1:nSubj
    cond2use_pha{subji} = cond2use_pha{subji}(:,:,mlim);
end

cfg_real                =   [];
cfg_real.pval           =   0.05; 
cfg_real.cond2use_pha   =   cond2use_pha;
cfg_real.rsa2cmp        =   perTrialRSA;


[out_real]         =   real_diff_pha_3D(cfg_real);


if runPerm
    cfg_perm                   =       [];
    cfg_perm.nPerm             =       nPerm;            
    cfg_perm.pval              =       0.05;
    cfg_perm.out_real          =       out_real;
    cfg_perm.cond2use_pha      =       cond2use_pha;
    cfg_perm.rsa2cmp           =       perTrialRSA;
    cfg_perm.savePerm          =       1;

    [out_perm rsa_pha_corr]    = myPerm_pha_3D (cfg_perm);
end

toc;

%% get rank for stats
max_clust_sum_ranked = out_perm.max_clust_sum;
max_clust_sum_ranked(:,2) =  tiedrank (out_perm.max_clust_sum);

[c index] = min(abs(max_clust_sum_ranked(:,1)- out_perm.max_clust_sum_real));
1000 - max_clust_sum_ranked(index, 2)
%x = find(max_clust_sum_ranked(:,2) == index)
%% plot grand average phase and rsa correlations
%times = 1:40;
times = 0:.1:1.7;
%times = 0:.1:3.4;
%min_freq1 = 1;max_freq1 = 29;num_frex1 = 29;min_freq2 = 30;max_freq2 = 100;num_frex2 = 15;baseline_window = [ -500 0]; frex1 = linspace(min_freq1,max_freq1,num_frex1);frex2 = linspace(min_freq2,max_freq2,num_frex2);frex = [frex1 frex2]; num_frex = size(frex, 2);
frex = 1:100;

cor2plot = out_real.rsa_pha_corrMean;
[h p ci t] = ttest (out_real.rsa_pha_corr, 0, 'alpha', 0.05); 
h = squeeze(h);
t = squeeze(t.tstat);
clustInfoReal = bwconncomp(h);
numPixReal = cellfun(@numel,clustInfoReal.PixelIdxList);
[max_clust_infoR maxiR] = max([ cellfun(@numel,clustInfoReal.PixelIdxList) ]); % the zero accounts for empty maps
clear all_clust_tsum_real;
for ci = 1:length(clustInfoReal.PixelIdxList)
    all_clust_tsum_real(ci) = sum (t(clustInfoReal.PixelIdxList{ci}));
end
all_clust_tsum_real = all_clust_tsum_real';
all_clust_tsum_real(:,2) =  numPixReal;


figure(1);
subplot(211)
contourf(times,frex, cor2plot, 40, 'linecolor', 'none');hold on
%contour (times, frex, h, 1, 'lineWidth', 1, 'linecolor', 'k');
%plot([0 0],get(gca,'ylim'),'k');axis square; hold on;
plot([0 0],get(gca,'ylim'),'k');axis square; hold on; 
%set(gca, 'xlim', [-.5 3], 'clim', [-0.25 0.25]); colorbar %, 'ylim', [1 100],
set(gca, 'clim', [-0.25 0.25]); colorbar %, 'ylim', [1 100],

subplot(212)
%contourf(times, frex, t, 40, 'linecolor', 'none');hold on
lab2plot = {'20', '40', '60', '80', '100'};
lab2plot = fliplr(lab2plot);
imagesc(times, frex, flipud(t));hold on
set (gca, 'ytick', 1:20:100, 'yticklabel', lab2plot);
h = zeros(100, size(h,2));
h(out_real.clustInfoReal.PixelIdxList{1}) = 1;
contour (times, frex, flipud(h), 1, 'lineWidth', 2, 'linecolor', 'k');
%contour (times, frex, h, 1, 'lineWidth', 2, 'linecolor', 'k');

%set(gca, 'xlim', [-.5 3], 'clim', [-3 3]); colorbar %, 'ylim', [1 100],
set(gca,  'clim', [-3 4], 'FontSize', 8); axis square; colorbar %, 'ylim', [1 100],
load myCmap;
colormap (myCmap);

export_fig(1, '_ispc_rsa_allTS_corr.png','-transparent', '-r300');


%% plot phase and rsa correlations for each subject
times = -3.75:.1:3.75;
subj = 6;
cor2plot = squeeze(rsa_pha_corr(subj,:,:));
contourf(cor2plot, 40, 'linecolor', 'none'); axis square;hold on
%contour (times, frex, h, 1, 'lineWidth', 1, 'linecolor', 'k');
%plot([0 0],get(gca,'ylim'),'k');hold on;
%set(gca, 'xlim', [-.5 3], 'clim', [-3 3]); colorbar %, 'ylim', [1 100],
%set(gca, 'xlim', [-.5 3], 'clim', [-.05 .05], 'ylim', [1 50]); %

%% remove sub2
SI_pha{2} = [];
SI_rsa_H{2} = []; SI_rsa{2} = [];
SI_pha = SI_pha(~cellfun('isempty',SI_pha));
SI_rsa_H = SI_rsa_H(~cellfun('isempty',SI_rsa_H));
SI_rsa = SI_rsa(~cellfun('isempty',SI_rsa));
SI_coh_i{2} = [];
SI_coh_i = SI_coh_i(~cellfun('isempty',SI_coh_i));
disp ('subj removed');

%% tROI phase sync analysis (1 value per trial)
clear, close all
tic

sublist = dir('*_coh_inROI.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    load(sublist{subji});
    % baseline normalization
    %all{subji,1} = ispcRoi; 
    all{subji,1} = ispcRoi - ispcBaseline; 
    %all{subji,1} = ispcRoi - mean(ispcBaseline);
    %all{subji,1} = bsxfun(@minus, ispcRoi, mean(ispcBaseline, 1));
   
end


filename = [sublist{subji}(4:end-8)];
eval([filename '= all;']);

% meanCOH_tmp     = cell2mat(cellfun(@mean, allM, 'un',0)); 
% meanCOH         = mean(meanCOH_tmp);
% stdCOH          = std(meanCOH_tmp, [], 1);
% seCOH           = stdCOH/sqrt(6);

cd ..
save (filename, filename);

toc


%% test
A = randi([0 100], 40, 100) / 100;
B = randi([0 100], 40, 100) / 100;
C = mean(B);
D = A-B;
E = A-C;
T1 = mean (D);
T2 = mean (E);

%% LOAD DATA and remove sub2
%clearvars;
load SI_coh_i; load perTrialRSA_H_SI;
SI_coh_i{2} = []; SI_coh_i = SI_coh_i(~cellfun('isempty',SI_coh_i));
perTrialRSA{2} = []; perTrialRSA = perTrialRSA(~cellfun('isempty',perTrialRSA));
disp ('subj removed');
%% calculate mean rsa in cluster

takeClust       = 0;
cond2use_pha    = SI_coh_i;
runperm         = 0;
nPerm           = 1000;

%calculate real differences
cfg_real                =   [];
cfg_real.pval           =   0.05; 
cfg_real.cond2use_pha   =   cond2use_pha;
cfg_real.rsa2cmp        =   perTrialRSA;


[out_real]         =   real_diff_2D(cfg_real);


if runperm
    cfg_perm                   =       [];
    cfg_perm.nPerm             =       nPerm;            
    cfg_perm.pval              =       0.05;
    cfg_perm.out_real          =       out_real;
    cfg_perm.cond2use_pha      =       cond2use_pha;
    cfg_perm.rsa2cmp           =       perTrialRSA;
    cfg_perm.savePerm          =       1;

    out_perm = myPerm2d (cfg_perm);
end


%%plot
figure(1); set(gcf,'Position', [0 0 900 700]); 
left_color = [0 0 0];
right_color = [0 0 0];
set(1,'defaultAxesColorOrder',[left_color; right_color]);
frex = 1:100;

meanCOH_tmp     = cell2mat(cellfun(@mean, cond2use_pha, 'un',0)); 
meanCOH         = mean(meanCOH_tmp);
stdCOH          = std(meanCOH_tmp, [], 1);
seCOH           = stdCOH/sqrt(6);

shadedErrorBar(frex, out_real.rsa_pha_corrMean, out_real.rsa_pha_var, 'r', 1); hold on;
plot (get(gca, 'xlim'), [0 0], 'k:');
h1 = out_real.sigMH_real'; h2p = h1';h2p(h2p==1) = .35;h2p(h2p==0) = NaN;
plot (h2p, 'k-', 'LineWidth', 2);

set(gca, 'ylim', [-.35 .4]);
axesHandles = findall(0, 'type', 'axes');
set(axesHandles, 'FontSize', 30); %, 'xticklabel',{'s1', 's2', 's3', 's4', 's5', 's6', 's7'}
xlabel(axesHandles(1),{'Frequency'})
ylabel(axesHandles(1),{['Correlation of phase synchronization' newline ...
    'and hippocampal reinstatement']})
yyaxis right
ylabel(axesHandles(1),{'Phase synchronization'}) 
set(gca, 'ylim', [-.1 .1]);
shadedErrorBar(frex,meanCOH , seCOH, 'k-', 1); hold on;


export_fig(1, '_ispc_rsa_tROI_corr.png','-transparent', '-r300');



%% get rank for stats
max_clust_sum_ranked = out_perm.max_clust_sum;
max_clust_sum_ranked(:,2) =  tiedrank (out_perm.max_clust_sum);

[c index] = min(abs(max_clust_sum_ranked(:,1)- (abs(out_perm.max_clust_sum_real))));
1000 - max_clust_sum_ranked(index, 2)
%x = find(max_clust_sum_ranked(:,2) == index)

%% plot 
%min_freq1 = 1;max_freq1 = 29;num_frex1 = 29;min_freq2 = 30;max_freq2 = 100;num_frex2 = 15;baseline_window = [ -500 0]; frex1 = linspace(min_freq1,max_freq1,num_frex1);frex2 = linspace(min_freq2,max_freq2,num_frex2);frex = [frex1 frex2]; num_frex = size(frex, 2);
frex = 1:100;

rsa_pha_corrMean = mean (rsa_pha_corr,1);
[h p ci t] = ttest (rsa_pha_corr);

figure(1);
subplot (121)
imagesc (flipud(rsa_pha_corr'));colorbar;
%imagesc (rsa_pha_corr');colorbar;
subplot (122)
im2plot = rsa_pha_corrMean';
im2plotf = flipud(im2plot );
imagesc (im2plotf); hold on;
h_copy = h;
h = h';h1 = h;
h1(:,2) = h; %h1 = flipud(h1);
contour(h1, 1, 'lineWidth', 3, 'linecolor', 'k');hold on
set(gca, 'xlim', [1.2 1.5], 'clim', [-.3 .15]);colorbar; 

frex4lab = fliplr(frex); frex4lab = frex4lab(1:4:end);
%frex4lab = frex; frex4lab = frex(1:4:end);
axesHandles = findall(0, 'type', 'axes');
set(axesHandles, 'ytick', 1:4:100, 'yticklabel',frex4lab,  'FontSize', 14, ...
    'xtick', 1:1:7 ); %, 'xticklabel',{'s1', 's2', 's3', 's4', 's5', 's6', 's7'}
ylabel(axesHandles(1),{'Frequency (Hz)'})
xlabel(axesHandles(1),{'Grand Average'})
ylabel(axesHandles(2),{'Frequency (Hz)'})
xlabel(axesHandles(2),{'Subject'})


export_fig(1, '_ispc_rsa_tROI_corr.png','-transparent', '-r300');

%% get rank (for stats)
max_clust_sum_ranked = out_perm.max_clust_sum;
max_clust_sum_ranked(:,2) =  tiedrank (out_perm.max_clust_sum);

[c index] = min(abs(max_clust_sum_ranked(:,1)- out_perm.max_clust_sum_real));
1000 - max_clust_sum_ranked(index, 2)
%x = find(max_clust_sum_ranked(:,2) == index)




%% LOAD DATA and exclude subjects from the saved contrasts
clearvars; 
load SI_H; load SI_LTC; load clustInfoReal_H; load clustInfoReal_LTC_HC_LC;
load SI_coh_i;
s = [2];
HC{s} = []; HC = HC(~cellfun('isempty',HC)); 
SI_H{s} = []; SI_H = SI_H(~cellfun('isempty',SI_H)); 
SI_LTC{s} = []; SI_LTC = SI_LTC(~cellfun('isempty',SI_LTC)); 
SI_coh_i{s} = []; SI_coh_i = SI_coh_i(~cellfun('isempty',SI_coh_i)); 
SISP_H{s} = []; SISP_H = SISP_H(~cellfun('isempty',SISP_H)); 
SISP_LTC{s} = []; SISP_LTC = SISP_LTC(~cellfun('isempty',SISP_LTC)); 
SIDR_H{s} = []; SIDR_H = SIDR_H(~cellfun('isempty',SIDR_H)); 
SIDR_LTC{s} = []; SIDR_LTC = SIDR_LTC(~cellfun('isempty',SIDR_LTC)); 



%% phase sync confidence analysis
%%import Phase sync
%can be used to import global and local RSA
%clear, close all
clearvars 
tic

sublist = dir('*_coh.mat'); 
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    load(sublist{subji});
    %chan2plot(subji)
    size(ispc)
    blne = [36:36]; %this blne works= [32:38];
    %ispcB = 100* bsxfun(@rdivide, ispc - mean(ispc(:,:,blne(1):blne(end)), 3), mean(ispc(:,:,blne(1):blne(end)), 3)); %pick first trial just to check
    ispcB = ispc - mean(ispc(:,:,blne(1):blne(end)), 3);
    all{subji,1} = ispcB; %when more than 1 electrode
    %all{subji,1} = rsaZ;
    
end

if exist ('timeBins')
    timeBins1 = timeBins (:, [1, end]);
end


cd .. % goes up one directory
filename = [sublist{subji}(4:end-8) '_pha'];
%filename = [sublist{subji}(5:end-12)];
eval([filename '= all;']);
save (filename, filename, '-v7.3');
%save (filename, filename);

toc

%% remove sub2
clearvars;
load HC_pha; load LC_pha
HC_pha{2} = [];
HC_pha = HC_pha(~cellfun('isempty',HC_pha));
disp ('subj removed');
%% plot mean for one condition 
all = LC_pha;
m_all = cell2mat(cellfun(@mean, all, 'un',0)); 
[h p ci t] = ttest (m_all, 0, 'alpha', 0.05); 
h = squeeze(h);
t = squeeze(t.tstat);
%44F
%min_freq1 = 1;max_freq1 = 29;num_frex1 = 29;min_freq2 = 30;max_freq2 = 100;num_frex2 = 15;baseline_window = [ -500 0]; frex1 = linspace(min_freq1,max_freq1,num_frex1);frex2 = linspace(min_freq2,max_freq2,num_frex2);frex = [frex1 frex2]; num_frex = size(frex, 2);
%100F
frex = 1:100;

grand_mean_all = squeeze(mean(m_all));
%times = -3.75:.1:3.75;
times = 1:76;

figure();
subplot (211);
contourf(times, frex, grand_mean_all, 40, 'linecolor', 'none'); axis square;hold on
%plot([0 0],get(gca,'ylim'),'k');hold on;
plot([41 41],get(gca,'ylim'),'k');hold on;
%set(gca, 'xlim', [-.5 3], 'clim', [-.025 .025]); colorbar; %colorbar, 'ylim', [1 100],
set(gca, 'xlim', [35 75], 'clim', [-.025 .025]); colorbar %, 'ylim', [1 100],

subplot(212)
contourf(times, frex, t, 40, 'linecolor', 'none'); axis square;hold on
contour (times, frex, h, 1, 'lineWidth', 1, 'linecolor', 'k');
%plot([0 0],get(gca,'ylim'),'k');hold on;
plot([41 41],get(gca,'ylim'),'k');hold on;
%set(gca, 'xlim', [-.5 3], 'clim', [-3 3]); colorbar %, 'ylim', [1 100],
set(gca, 'xlim', [35 75], 'clim', [-3 3]); colorbar %, 'ylim', [1 100],

export_fig(1, '_phase_sync.png','-transparent', '-r300');

%% take only first bin (figure 6A)
clearvars;
load HC_pha; load LC_pha; 

for subji = 1:6
    HC_bin(subji,:) = mean (squeeze( mean ( HC_pha{subji}(:,70:83,41), 2) ));
    LC_bin(subji,:) = mean (squeeze( mean ( LC_pha{subji}(:,70:83,41), 2) ));
end

[h p ci t] = ttest(HC_bin(1:5), LC_bin(1:5))

%%
clearvars;
load HC_pha; load LC_pha; 
for subji = 1:6
    HC_bin{subji} = squeeze( mean ( HC_pha{subji}(:,70:83,41) ));
    LC_bin{subji} = squeeze( mean ( LC_pha{subji}(:,70:83,41) ));
end

%% HC LC Phase sync 
%2 bar
data.data = [HC_bin LC_bin];

figure(1); set(gcf,'Position', [0 0 500 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(hb, 'linestyle','none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LC  '}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.075 .1] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,1), data.data(:,2));

disp (['p = ' num2str(p) '   t= ' num2str(t.tstat)]);

export_fig(1, '_jackknife_2bars.png','-transparent', '-r300');
close all;





%% Take until 1.7s
clearvars;
load HC_pha; load LC_pha; load clustInfoReal_LTC_SI;
HC_pha{2} = []; HC_pha = HC_pha(~cellfun('isempty',HC_pha));


all_H = HC_pha; 
all_L = LC_pha;
m_all_H = cell2mat(cellfun(@mean, all_H, 'un',0)); 
m_all_L = cell2mat(cellfun(@mean, all_L, 'un',0)); 

frex = 1:100;
mlim = 41:58; % until the end of the LTC conf cluster
% apply limits 
nSubj = length(all_H);
for subji = 1:nSubj
    m_all_H_lim(subji, :, :) = m_all_H(subji,:,mlim);
    m_all_L_lim(subji, :, :) = m_all_L(subji,:,mlim);
end

grand_mean_all_H = squeeze(mean(m_all_H_lim));
grand_mean_all_L = squeeze(mean(m_all_L_lim));

for subji = 1:nSubj
    
    m_perS_H(subji, :) = mean (m_all_H_lim(subji, clustInfoReal.PixelIdxList{1}) ) ;
    m_perS_L(subji, :) = mean (m_all_L_lim(subji, clustInfoReal.PixelIdxList{1}) ) ;
    
end
 

test = zeros(size(grand_mean_all_H,1), size(grand_mean_all_H, 2)); 
test(clustInfoReal.PixelIdxList{1}) = 1;

meanInClust_H = mean (grand_mean_all_H(clustInfoReal.PixelIdxList{1}) ); 
meanInClust_L = mean (grand_mean_all_L(clustInfoReal.PixelIdxList{1}) ); 

%% HC _ LC phase sync
%2 bar

data.data = [m_perS_H m_perS_L];

figure(1); set(gcf,'Position', [0 0 500 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(hb, 'linestyle','none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LC  '}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.15 .1] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['p = ' num2str(p) '   t= ' num2str(t.tstat)]);

export_fig(1, '_HC_LC_phaSync_2bars.png','-transparent', '-r300');
close all;



%%
figure();
subplot (211);
contourf(mlim, frex, grand_mean_all, 40, 'linecolor', 'none'); axis square;hold on
%plot([0 0],get(gca,'ylim'),'k');hold on;
%plot([41 41],get(gca,'ylim'),'k');hold on;
%set(gca, 'xlim', [-.5 3], 'clim', [-.025 .025]); colorbar; %colorbar, 'ylim', [1 100],
%set(gca, 'xlim', [35 75], 'clim', [-.025 .025]); colorbar %, 'ylim', [1 100],

subplot(212)
contourf(test)
%contourf(mlim, frex, t, 40, 'linecolor', 'none'); axis square;hold on
%contour (mlim, frex, h, 1, 'lineWidth', 1, 'linecolor', 'k');



export_fig(1, '_phase_sync.png','-transparent', '-r300');



%% multiple regression analysis
clearvars;
load multipleRegressionAnalysis;

% create table for each subject
for subji = 1:6
    ERS_H               =   perTrialRSA_H{subji}'; 
    ERS_LTC             =   perTrialRSA_LTC{subji}';
    PLV                 =   mean (SI_coh_i{subji}(:,70:83),2);
    tbl{subji}          =   table (ERS_H, ERS_LTC, PLV,...
                            'VariableNames',{'ERS_H','ERS_LTC','PLV'});
%     mdl{subji}          = fitlm(tbl{subji},'ResponseVar','ERS_LTC',...
%                             'PredictorVars',{'ERS_H'}); %'RobustOpts','on'
    mdl{subji}          = fitlm(tbl{subji},'interactions','ResponseVar','ERS_LTC',...
                            'PredictorVars',{'ERS_H','PLV'}); %'RobustOpts','on'
    mdl{subji}
    beta_ERS_H(subji)   = mdl{subji}.Coefficients.Estimate(2);
    beta_PLV(subji,:)     = mdl{subji}.Coefficients.Estimate(3);
    beta_Inter(subji)   = mdl{subji}.Coefficients.Estimate(4);
    Rsquared(subji)     = mdl{subji}.Rsquared.Ordinary;
end


%% multiple regression analysis 2
clearvars;
load multipleRegressionAnalysis;

% create table for each subject
for subji = 1:6
    ERS_H               =   perTrialRSA_H{subji}'; 
    ERS_LTC             =   perTrialRSA_LTC{subji}';
    PLV                 =   mean (SI_coh_i{subji}(:,70:83),2);
    X                   =   [ones(size(ERS_H)) ERS_H PLV ERS_H.*PLV];
    %X                   =   [ones(size(ERS_H)) ERS_H PLV];
    %X                   =   [ones(size(ERS_H)) PLV];
    y                   =   ERS_LTC;
    b                   =   regress(y,X);
    a                   =   X\y;
    Y                   =   X*a;
    MaxErr(subji)       = max(abs(Y - y));
end


MaxErr = MaxErr';
mean(MaxErr)

%%
scatter3(ERS_H,PLV,y,'filled')
hold on
x1fit = min(ERS_H):0.01:max(ERS_H);
x2fit = min(PLV):0.01:max(PLV);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('HC-RSA')
ylabel('PLV')
zlabel('LTC-RSA')
view(50,10)
hold off

%% multiple regression analysis 3
%clearvars;
load multipleRegressionAnalysis;

% create table for each subject
for subji = 1:6
    ERS_H                   =   perTrialRSA_H{subji}'; 
    ERS_LTC                 =   perTrialRSA_LTC{subji}';
    PLV                     =   mean (SI_coh_i{subji}(:,70:83),2);
    s{subji}                =   regstats(ERS_LTC, [PLV ERS_H]);
    PLV_factor(subji, :)    =   s{subji}.beta(2);
    H_ERS_factor(subji, :)  =   s{subji}.beta(3);
    beta{subji}             =   s{subji}.beta;
    rSquared{subji}         =   s{subji}.rsquare;
end



%%
scatter3(ERS_H,PLV,y,'filled')
hold on
x1fit = min(ERS_H):0.01:max(ERS_H);
x2fit = min(PLV):0.01:max(PLV);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('HC-RSA')
ylabel('PLV')
zlabel('LTC-RSA')
view(50,10)
hold off

%% partial correlation
clearvars;
load multipleRegressionAnalysis;

for subji = 1:6
    ERS_H                           =   perTrialRSA_H{subji}'; 
    ERS_LTC                         =   perTrialRSA_LTC{subji}';
    PLV                             =   mean (SI_coh_i{subji}(:,70:83),2);
    [r_corr1(subji),p_corr1(subji)] =   corr(PLV, ERS_LTC, 'Type','s'); %Spearman
    [r_corr2(subji),p_corr2(subji)] =   corr(ERS_H, PLV, 'Type','s'); %Spearman
    [r_part1(subji),p_part1(subji)] =   partialcorr(ERS_H, ERS_LTC, PLV, 'Type','s'); %Spearman
    [r_part2(subji),p_part2(subji)] =   partialcorr(PLV, ERS_LTC, ERS_H, 'Type','s'); %Spearman

end

[h p ci t] = ttest(r_corr2)
[h p ci t] = ttest(r_part2)

%%
clearvars
load carsmall
x1 = Weight;
x2 = Horsepower;    % Contains NaN data
y = MPG;

X = [ones(size(x1)) x1 x2 x1.*x2];
b = regress(y,X)    % Removes NaN data


%%
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):100:max(x1);
x2fit = min(x2):10:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Weight')
ylabel('Horsepower')
zlabel('MPG')
view(50,10)
hold off



%%
clearvars;
x1 = [.2 .5 .6 .8 1.0 1.1]';
x2 = [.1 .3 .4 .9 1.1 1.4]';
y  = [.17 .26 .28 .23 .27 .24]';

X = [ones(size(x1))  x1  x2];

a = X\y



%%
data.data = [beta_PLV'  beta_ERS_H' beta_Inter'];

figure(1); set(gcf,'Position', [0 0 500 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2 3], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(hb, 'linestyle','none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3],'XTickLabel',{'PLV  ', 'H_ERS  ', 'Interaction  '}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 4], 'ylim', [-.085 7.75] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,3));

disp (['t= ' num2str(t.tstat) '   p = ' num2str(p)   ]);

export_fig(1, '_regressors.png','-transparent', '-r300');
close all;


%% COORDINATED REINSTATEMENT IN SPECIFIC CLUSTERS BY PARTIALLING OUT ISPC
%clearvars;
load clustInfoReal_H ; load clustInfoReal_LTC_HC_LC;

saveimg = 1;
takeClusters = 0;
clear all_trials_RSA_H all_trials_RSA_LTC;
all_trials_RSA_H    = SI_H;
all_trials_RSA_LTC  = SI_LTC;
all_trials_ISPC     = SI_coh_i;
freqRange           = 70:82;

%40
pixel_LTC = 12;% HC-LC_in_LTC = 12;pixel = 17; LTC_SIDI = 17-20
mlim = 36:75;

%%time matched 
pixel_H = 6; %H_cong_inc = 6; %1 for matched at encoding %
roiX_H  = 25:31; roiY_H  = 6:11; % 
roiX_LTC = 25:31; roiY_LTC  = 13:23; % 

%%most significant 
% pixel_H = 14; %use biggestSmallest H cluster. H_cong_inc = 6; %1 for matched at encoding %
% roiX_H  = 23:28; roiY_H  = 7:11; % 
% roiX_LTC = 26:31; roiY_LTC  = 14:22; % 

clear perTrialRSA_H perTrialRSA_LTC subjt_H subjt_LTC;
for subji = 1:length(all_trials_RSA_H)
    allTRSA_H = all_trials_RSA_H{subji}; 
    allTRSA_H = allTRSA_H (:, mlim, mlim);
    allTRSA_LTC = all_trials_RSA_LTC{subji};
    allTRSA_LTC = allTRSA_LTC (:, mlim, mlim);
    for ti = 1:size(allTRSA_H, 1)
        tMat = allTRSA_H(ti, roiX_H, roiY_H); %square
        if (takeClusters)
            tMat = allTRSA_H(ti, clustInfoReal_H.PixelIdxList{pixel_H}); %cluster
        end
        t = tMat(:);   t = mean (t); 
        subjt_H(ti,:) = t;
        perTrialRSA_H{subji} = subjt_H'; 
    end
    for ti = 1:size(allTRSA_LTC, 1)
        t = allTRSA_LTC(ti, roiX_LTC, roiY_LTC);  %square
        if (takeClusters)
            t = allTRSA_LTC(ti, clustInfoReal_LTC.PixelIdxList{pixel_LTC}); %cluster
        end
        t = t(:);   t = mean (t); %cluster
        subjt_LTC(ti,:) = t;
        perTrialRSA_LTC{subji} = subjt_LTC'; 
    end
    
    clear subjt_H subjt_LTC;
end
perTrialRSA_H = perTrialRSA_H'; perTrialRSA_LTC = perTrialRSA_LTC';


% define frequency range for ispc
clear perTrial_ISPC; 
for subji = 1:length(all_trials_ISPC)
    perTrial_ISPC{subji} = squeeze( mean (all_trials_ISPC{subji}(:,freqRange),2 )) ;
end
perTrial_ISPC = perTrial_ISPC';


%%CORRELATE RSA per trial in H and LTC
clear rsa;
for subji = 1:length(perTrialRSA_H)
    x = perTrialRSA_H{subji}(:);
    y = perTrialRSA_LTC{subji}(:);
    z = perTrial_ISPC{subji}(:);
    %[r,p] = corr(x, y, 'Type','s'); 
    [r,p] = partialcorr(x, y, z, 'Type','s'); %Spearman
    rsa(subji,1) = r; rsa(subji,2) = p;
    meanX(subji) = mean(x);
    meanY(subji) = mean(y);
end

rsaZ = atanh(squeeze(rsa));
rsaZ(:,1)'
% t-test
[h p ci t] = ttest (rsaZ(:,1)); %1 means Rhos 
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

clear test;
test = squeeze(mean(allTRSA_LTC, 1));
test(1:end) = 0.2;
%squares
test(roiX_H,roiY_H,:) = .55;
test(roiX_LTC,roiY_LTC,:) = .55;

%pixels
test(clustInfoReal_H.PixelIdxList{pixel_H}) = 1; test(clustInfoReal_H.PixelIdxList{pixel_H}) + 1;
test(clustInfoReal_LTC.PixelIdxList{pixel_LTC}) = 1;
test = flipud(test); 
figure(1);imagesc(test);axis square;hold on;
%title(['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
%remove labels
axesHandles = findall(0, 'type', 'axes');
for i=1:length(axesHandles)
   set (axesHandles(i), 'visible', 'off'); % remove labels 
end



if saveimg
titFig = 'current_fig';
export_fig (titFig, '-transparent', '-r200');
close all;
end


%% plot  bars
clear data;
data.data = out_real.rsa_pha_corr';

figure(1); set(gcf,'Position', [0 0 400 400]); 
mean_S = mean(data.data, 2);
hb = plot ([1:100], data.data, '.'); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',24);hold on;

h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1:100], 'FontSize', ...
    10, 'linew',2, 'xlim', [0 101], 'ylim', [-1 1] );
% set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'All  ', '1-15  ', '15-30  ',...
%     '30-100  ', '1-8  ', '9-29  '}, 'FontSize', ...
%     20, 'linew',2, 'xlim', [0 7], 'ylim', [-.1 .2] );
xticklabel_rotate;
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

export_fig(1, '_jackknife.png','-transparent', '-r300');
%close all;



%% remove sub 4 partialling out analysis
load SI_coh_H_ROI; load SI_coh_LTC_ROI; load SI_rsa_H; load SI_rsa_LTC;
subj2exc = 2;
SI_coh_H_ROI{subj2exc} = []; SI_coh_LTC_ROI{subj2exc} = []; 
SI_rsa_H{subj2exc} = []; SI_rsa_LTC{subj2exc} = [];
SI_coh_H_ROI = SI_coh_H_ROI(~cellfun('isempty',SI_coh_H_ROI));
SI_rsa_H = SI_rsa_H(~cellfun('isempty',SI_rsa_H));
SI_coh_LTC_ROI = SI_coh_LTC_ROI(~cellfun('isempty',SI_coh_LTC_ROI));
SI_rsa_LTC = SI_rsa_LTC(~cellfun('isempty',SI_rsa_LTC));
disp ('subj removed');

%% COORDINATED REINSTATEMENT BY PARTIALING OUT ISPC
%

saveimg = 1;
takeClusters = 0;
clear all_trials_RSA_H all_trials_RSA_LTC;
all_trials_RSA_H = SI_rsa_H;
all_trials_RSA_LTC = SI_rsa_LTC;
all_trials_phase = SI_coh_H_ROI; %> DEFINE HERE what retrieval time to use for ISPC

%40
pixel_LTC = 12;% HC-LC_in_LTC = 12;pixel = 17; LTC_SIDI = 17-20
mlim = 36:75;

%%time matched 
pixel_H = 6; %H_cong_inc = 6; %1 for matched at encoding %
roiX_H  = 25:31; roiY_H  = 6:11; % 
roiX_LTC = 25:31; roiY_LTC  = 14:22; % 


clear perTrialRSA_H perTrialRSA_LTC subjt_H subjt_LTC;
for subji = 1:length(all_trials_RSA_H)
    allTRSA_H = all_trials_RSA_H{subji}; 
    allTRSA_H = allTRSA_H (:, mlim, mlim);
    allTRSA_LTC = all_trials_RSA_LTC{subji};
    allTRSA_LTC = allTRSA_LTC (:, mlim, mlim);
    for ti = 1:size(allTRSA_H, 1)
        tMat = allTRSA_H(ti, roiX_H, roiY_H); %square
        if (takeClusters)
            tMat = allTRSA_H(ti, clustInfoReal_H.PixelIdxList{pixel_H}); %cluster
        end
        t = tMat(:);   t = mean (t); 
        subjt_H(ti,:) = t;
        perTrialRSA_H{subji} = subjt_H'; 
    end
    for ti = 1:size(allTRSA_LTC, 1)
        t = allTRSA_LTC(ti, roiX_LTC, roiY_LTC);  %square
        if (takeClusters)
            t = allTRSA_LTC(ti, clustInfoReal_LTC.PixelIdxList{pixel_LTC}); %cluster
        end
        t = t(:);   t = mean (t); %cluster
        subjt_LTC(ti,:) = t;
        perTrialRSA_LTC{subji} = subjt_LTC'; 
    end
    
    clear subjt_H subjt_LTC;
end


perTrialRSA_H = perTrialRSA_H'; perTrialRSA_LTC = perTrialRSA_LTC';

%%CORRELATE RSA per trial in H and LTC by PARTIALING OUT ISPC
clear rsa rsaPO;
for freqi = 1:44
    freqi
    for subji = 1:length(perTrialRSA_H)
        x = perTrialRSA_H{subji}(:);
        y = perTrialRSA_LTC{subji}(:);
        z = all_trials_phase{subji}(:,freqi);
        [r p] = corr(x, y, 'Type','s'); %Spearman
        
        [rp pp] = partialcorr(x, y, z, 'Type','s'); %Spearman
        rsaPO(subji,freqi) = r - rp; %rsaPO(subji,freqi, 2) = p;
        %rsaPO(subji,freqi) = rp; %rsaPO(subji,freqi, 2) = p;
        
    end
end

%% plot  bars for partial out analysis
clear data;
data.data = rsaPO';

figure(1); set(gcf,'Position', [0 0 400 400]); 
mean_S = mean(data.data, 2);
hb = plot ([1:44], data.data, '.'); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',24);hold on;

h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1:44], 'FontSize', ...
    10, 'linew',2, 'xlim', [0 45], 'ylim', [-.1 .2] );
% set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'All  ', '1-15  ', '15-30  ',...
%     '30-100  ', '1-8  ', '9-29  '}, 'FontSize', ...
%     20, 'linew',2, 'xlim', [0 7], 'ylim', [-.1 .2] );
xticklabel_rotate;
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

export_fig(1, '_jackknife.png','-transparent', '-r300');
%close all;



%% plot 
min_freq1 = 1;max_freq1 = 29;num_frex1 = 29;min_freq2 = 30;max_freq2 = 100;num_frex2 = 15;baseline_window = [ -500 0]; frex1 = linspace(min_freq1,max_freq1,num_frex1);frex2 = linspace(min_freq2,max_freq2,num_frex2);frex = [frex1 frex2]; num_frex = size(frex, 2);

rsaPOMean = mean (rsaPO,1);

figure(1);
subplot (121)
imagesc (flipud(rsaPO'));colorbar;
%imagesc (rsa_pha_corr');colorbar;
subplot (122)
im2plot = rsaPOMean';
im2plotf = flipud(im2plot );
imagesc (im2plotf); hold on;

frex4lab = fliplr(frex); frex4lab = frex4lab(1:4:end);
%frex4lab = frex; frex4lab = frex(1:4:end);
axesHandles = findall(0, 'type', 'axes');
set(axesHandles, 'ytick', 1:4:44, 'yticklabel',frex4lab,  'FontSize', 14, ...
    'xtick', 1:1:7 ); %, 'xticklabel',{'s1', 's2', 's3', 's4', 's5', 's6', 's7'}
ylabel(axesHandles(1),{'Frequency (Hz)'})
xlabel(axesHandles(1),{'Grand Average'})
ylabel(axesHandles(2),{'Frequency (Hz)'})
xlabel(axesHandles(2),{'Subject'})


export_fig(1, '_rsa_pOUT_ISPC.png','-transparent', '-r300');

%%
rsaZ = atanh(squeeze(rsa));
rsaZ(:,1)'
% t-test
[h p ci t] = ttest (rsaZ(:,1)); %1 means Rhos 
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

clear test;
test = squeeze(mean(allTRSA_LTC, 1));
test(1:end) = 0.2;
%squares
test(roiX_H,roiY_H,:) = .55;
test(roiX_LTC,roiY_LTC,:) = .55;

%pixels
test(clustInfoReal_H.PixelIdxList{pixel_H}) = 1; test(clustInfoReal_H.PixelIdxList{pixel_H}) + 1;
test(clustInfoReal_LTC.PixelIdxList{pixel_LTC}) = 1;
test = flipud(test); 
figure(1);imagesc(test);axis square;hold on;
%title(['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
%remove labels
axesHandles = findall(0, 'type', 'axes');
for i=1:length(axesHandles)
   set (axesHandles(i), 'visible', 'off'); % remove labels 
end



if saveimg
titFig = 'current_fig';
export_fig (titFig, '-transparent', '-r200');
close all;
end



%% JACKKNIFE IN CLUSTER ANALYSIS
% first produce the *a. files (mean in cluster code section (save m_roi_cond1) 
%% MEAN IN CLUSTER OF INTEREST
clustInfoReal = clustInfoReal_LTC;

rsa_cond1 = cell2mat(cellfun(@mean, SI, 'un',0)); 

takeclust = 1;
%roiX  = 25:31; roiY  = 6:11; % % Hippocampus
roiX = 25:31; roiY = 14:22; % LTC
pixel = 12;%H_cong_inc = 6; % LTC_HC_LC = 12 % LTC_SIDI = 16-19

mlim = 36:75;
rsa_cond1 = rsa_cond1(:,mlim,mlim);

nSubj = size(rsa_cond1, 1);
if takeclust
    ROI_cond1 = rsa_cond1(:,clustInfoReal.PixelIdxList{pixel}); 
    m_roi_cond1 = mean(ROI_cond1, 2);
else
    ROI_cond1 = rsa_cond1(:,roiX, roiY); 
    ROI_cond1 = reshape (ROI_cond1, [nSubj size(roiX,2) * size(roiY, 2)]);
    m_roi_cond1 = mean(ROI_cond1, 2);
end

test = squeeze(mean(rsa_cond1, 1));
test(1:end) = 0;
test(roiX,roiY,:) = .5;
test(clustInfoReal.PixelIdxList{pixel}) = 1;
test = flipud(test);
figure(1);imagesc(test);axis square;hold on;
%figure(1);contourf(test, 40, 'linecolor', 'none');axis square;hold on;
% plot([mfL+0.5 mfL+0.5],get(gca,'ylim'),'w', 'LineWidth', 2); 
% plot(get(gca,'xlim'), [bins - (mfL-0.5) bins - (mfL-0.5)],'w', 'LineWidth', 2);
% plot(get(gca,'xlim'), [bins+.5 .5],'w', 'LineWidth', 2); % diagonal
% set (gca, 'clim', [0 1]);
%remove labels
axesHandles = findall(0, 'type', 'axes');
for i=1:length(axesHandles)
   set (axesHandles(i), 'visible', 'off'); % remove labels 
end


[h p ci t] = ttest (m_roi_cond1)
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

% set (gcf, 'Position', [200 200 400 450]);
% export_fig(2, 'plotting.png','-transparent', '-r300');
%close all;   





%% 
clearvars 
tic

sublist = dir('*a.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    load(sublist{subji});
    %chan2plot(subji)
    
    
    all(subji,:) = m_roi_cond1; %when more than 1 electrode
    %all{subji,1} = rsaZ;
    
end

toc
sublist = sublist';


%% exclude subject 2
all (:,2) = [];

%% plot  4 bar 
clear testminall data;
close all;
%testminall = (all(1,:) -all) %/ all(1,:); %testminall(:,) = [];
%data.data = testminall';
for subji = 1:6
    for freqi = 1:5
        all2save(subji) = (all(1, subji));
        testminall(freqi, subji) = ( (all(1, subji) - all(freqi, subji) )  ) ;
        %testminall(freqi, subji) = ( (all(1, subji) - all(freqi, subji) )  / all(1, subji) );      
    end 
end


data.data = testminall;


figure(1); set(gcf,'Position', [0 0 400 400]); 
mean_S = mean(data.data, 2);
hb = plot ([1 2 3 4 5], data.data, '.'); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',24);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'  ', '1-9 Hz  ',...
    '10-19 Hz  ', '20-29 Hz  ', '30-100Hz  '}, 'FontSize', ...
    18, 'linew',2, 'xlim',[1.1 5.9], 'ylim',[-0.075 .075] );%[-2.5 2.5]
xticklabel_rotate;
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

export_fig(1, '_jackknife.png','-transparent', '-r300');
%close all;



%% compare 1-8 in H and LTC
close all; 
freq = 2; %2 = theta; 3 = alpha ..

H = data_H.data(freq,[1 3:7]); H = H';
LTC = data_LTC.data(freq,[1 3:7]); LTC = LTC';
data.data = [H LTC];

figure(1); set(gcf,'Position', [0 0 500 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LTC  '}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.02 .075] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,1));

disp (['p = ' num2str(p) '   t= ' num2str(t.tstat)]);

export_fig(1, '_jackknife_theta.png','-transparent', '-r300');
close all;

%% compare difference all - (all - 1-100Hz) in H and LTC
%2 bar
%clear all; close all; 
load _m_roi_diff_H; load _m_roi_diff_LTC_HC_LC;

%data.data = [m_roi_diff_H([1 3:7],:) m_roi_diff_LTC_SIDI_C2([1 3:7],:)];
%data.data = [m_roi_diff_H([1 3:7],:) m_roi_diff_LTC_HC_LC];
%data.data = [m_roi_diff_all_H([1 3:7],:) m_roi_diff_all_LTC];
data.data = [m_roi_diff_all_minus_H([1 3:7],:) m_roi_diff_all_minus_LTC];

figure(1); set(gcf,'Position', [0 0 500 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(hb, 'linestyle','none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LTC  '}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.085 .15] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,2));

disp (['p = ' num2str(p) '   t= ' num2str(t.tstat)]);

export_fig(1, '_jackknife_2bars.png','-transparent', '-r300');
close all;


%% compare difference all - (all - 1-100Hz) in H and LTC
%4 bars %remove _r to get original diff plot
clear all; close all; 
load _m_roi_diff_H; load _m_roi_diff_LTC_HC_LC;
load _m_roi_diff_LTC_SI_DI_C1; load _m_roi_diff_LTC_SI_DI_C2;

data.data = [m_roi_diff_H([1 3:7],:) m_roi_diff_LTC_SIDI_C1([1 3:7],:)  ...
             m_roi_diff_LTC_SIDI_C2([1 3:7],:) m_roi_diff_LTC_HC_LC];

figure(1); set(gcf,'Position', [0 0 600 300]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2 3 4], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(hb, 'linestyle','none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'HC  ', 'HC  ', 
    'C1  ', 'C2 '}, ...
    'FontSize', 15, 'linew',2, 'xlim', [0.25 4.75], 'ylim', [-.04 .08] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
ylabel({['Encoding-retrieval' newline 'similarity reduction (rho)']})
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,4));

disp (['p = ' num2str(p) '   t= ' num2str(t.tstat)]);

export_fig(1, '_jackknife_4bars_HCSIDI.png','-transparent', '-r300');
close all;



%% compare difference all - (all - 1-100Hz) in H and LTC
% in 4 bars
clear all; close all; 
load m_roi_diff_H_SISP-SIDR_1-100Hz; load m_roi_diff_H_SISP-SIDR_10-100Hz;
load m_roi_diff_LTC_HC-LC_1-100Hz; load m_roi_diff_LTC_HC-LC_10-100Hz;

data.data = [m_roi_diff_all_H([1 3:7],:) m_roi_diff_all_minus_H([1 3:7],:),...
             m_roi_diff_all_LTC, m_roi_diff_all_minus_LTC];

figure(1); set(gcf,'Position', [0 0 500 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2 3 4], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb, 'linestyle','none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4],'XTickLabel', ...
    {'1-100Hz  ', '10-100Hz  ', '1-100Hz  ', '10-100Hz  '}, ...
    'FontSize', 20, 'linew',2, 'xlim', [0 5], 'ylim', [-.085 .175] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,3), data.data(:,4));

disp (['p = ' num2str(p) '   t= ' num2str(t.tstat)]);

export_fig(1, '_jackknife_4bars_bothCond.png','-transparent', '-r300');
close all;


%% anova
data4anova = data.data;
data4an = reshape (data4anova, [12 2]);

[p,tbl] = anova2(data4an,6);

%% anova repeated measures
%example
load fisheriris
t = table(species,meas(:,1),meas(:,2),meas(:,3),meas(:,4),...
'VariableNames',{'species','meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});

%% test
subjs = {'s1';'s3';'s4';'s5';'s6';'s7';};
t = table(subjs,data4anova(:,1),data4anova(:,2),data4anova(:,3),data4anova(:,4),...
'VariableNames',{'Subjects','meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});

rm = fitrm(t,'meas1-meas4~Subjects','WithinDesign',Meas);


%% plot  4 bar 
clear testminall data;
close all;

data.data = [all_H(1, :)' all_H(2,:)' all_LTC(1, :)' all_LTC(2,:)' ];


figure(1); set(gcf,'Position', [0 0 400 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2 3 4], data.data, '.'); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',24);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'1-100Hz  ', '10-100Hz  ',...
    '1-100 Hz  ', '10-100 Hz  '}, 'FontSize', ...
    18, 'linew',2, 'xlim',[0 5], 'ylim',[-0.05 .1] );%[-2.5 2.5]
xticklabel_rotate;
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

export_fig(1, '_jackknife.png','-transparent', '-r300');
%close all;


%% anova 
data4anova = data.data; data4anova1 = data.data;
data4an = reshape(data4anova, [12, 2]);

[p,tbl] = anova2(data4an,6);


%%
testminall = all - all(1,:);
figure();
imagesc(testminall );colorbar;
set (gca, 'clim', [-.03 .03]);

%% stats 

cond = data.data(:,6);

[h p ci t] = ttest (cond);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);




%% plot spectrograms

subj = 's1';
min_freq1 = 1;max_freq1 = 29;num_frex1 = 29; min_freq2 = 30;max_freq2 = 100;num_frex2 = 15;baseline_window = [ -800 -300]; frex1 = linspace(min_freq1,max_freq1,num_frex1); frex2 = linspace(min_freq2,max_freq2,num_frex2); frex = [frex1 frex2]; num_frex = size(frex, 2); 
tftimes = 1:2251;

ch = 1;
for ti = 1:size (all2all,1)
    for eri = 1:2
        figure(1);
        if (eri == 1) e = 'e'; else e = 'r';end
        contourf(tftimes,frex,squeeze(all2all(ti,eri, ch,:,:)),40,'linecolor','none');hold on;
        set(gca,'clim',[-15 15],'ydir','normal','xlim',xlim); 
        title ([num2str(ti) '_' e]);%colorbar;
        saveas (1,[num2str(ti) '_' e '_' subj 'Z_all.jpg']);
        close all;
    end
end



%% PLOT TRACES
clearvars 
tic

sublist = dir('*_traces.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    load(sublist{subji});
    %chan2plot(subji)
    all_traces{subji,1} =   oneListTraces_c; 
    all_power{subji,1} =    oneListPow; 
    
end



%% load all2all
tic     
clearvars -except ALLEEG
% Get a list of data files ready to be analyzed
sublist = dir('*all2all.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subno=1:length(sublist)
    disp (['File nº ' num2str(subno)]);
    load(sublist{subno});
    s_all{subno} = all2all; 
    size(all2all)
    if exist ('combs')
        combs_all{subno} = combs;
    end
    if exist ('chanNames')
        chanNames_all{subno} = chanNames;
    end
end
toc
%% CREATE BATCHES of N trials from all2all
Nsubj = length(s_all);
tic
batch_bin = 200;
for subji = 1:Nsubj
   all2all_b =  s_all{subji};
   chanNames = chanNames_all{subji};
   nTrials = size (all2all_b, 1)
   nBatches = ceil(nTrials/batch_bin)
   clear checkBatches;
   for batchi = 1:nBatches
       if batchi == 1 
        checkBatches{batchi} = batchi:batchi*batch_bin;
        all2all = all2all_b(batchi:batchi*batch_bin, :, :, :, :);
        filename = [sublist{subji} '_batch_'  num2str(batchi) '_all2all.mat'];
        save (filename, 'all2all', 'chanNames'); 
       else
           if (batchi) *batch_bin < nTrials
            checkBatches{batchi} = ((batchi-1)*batch_bin)+1:(batchi*batch_bin);
            all2all = all2all_b(((batchi-1)*batch_bin)+1:(batchi*batch_bin), :, :, :, :);
            filename = [sublist{subji} '_batch_'  num2str(batchi) '_all2all.mat'];
            save (filename, 'all2all', 'chanNames'); 
            else 
            disp ('final');
            checkBatches{batchi} = ((batchi-1)*batch_bin)+1:nTrials;
            all2all = all2all_b(((batchi-1)*batch_bin)+1:nTrials, :, :, :, :);
            filename = [sublist{subji} '_batch_'  num2str(batchi) '_all2all.mat'];
            save (filename, 'all2all', 'chanNames'); 
           end
       end
   end
   checkBatches
end
toc


%% LOAD BATCH FILES RECURSIVELY AND CREATE A FOLDER FOR EACH SUBJECT
% run in folder with all DI files 
% first create the folders based on the subjects's names
clear, close all
sublist = dir('*_rsa.mat'); sublist = {sublist.name};
str_tmp = 's0';
for filei=1:length(sublist)
   str = sublist{filei};
   if ~strcmp(str(1:3), str_tmp(1:3))
      mkdir(str(1:3))
      movefile(str, str(1:3))
   end
   if strcmp(str(1:3), str_tmp(1:3))
      movefile(str, str(1:3))
   end
   str_tmp = sublist{filei};
end



%% then go recursively through all of them to form the 
tic
folders = dir(); dirs = find(vertcat(folders.isdir));
folders = folders(dirs);

for foldi = 1:length(folders)
    direct = folders(foldi);
    clear allFiles;
    if strcmp(direct.name(1), 's')
        cd (direct.name)
        sublist = dir('*_rsa.mat');
        sublist = {sublist.name};
        disp (['number of batches -> ' num2str(length(sublist))]);
        for filei=1:length(sublist)
            load(sublist{filei});
            allFiles{filei} = rsaZ;    
        end

        rsaZ = cat (2, allFiles{:});

        cd .. % goes up one directory
        fname = strsplit(sublist{filei}, '_');
        filename = [fname{1} '_' fname{2} '_rsa'];
        save (filename, 'rsaZ');
    end
end


toc

%% LOAD BATCH FILES
clear, close all
tic
sublist = dir('*_rsa.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for filei=1:length(sublist)
    load(sublist{filei});
    all{filei} = rsaZ;    
end

rsaZ = cat (2, all{:});

cd .. % goes up one directory
filename = [sublist{filei}(1:5) '_rsa'];
save (filename, 'rsaZ');


toc

%% LOAD BATCH FILES RECURSIVELY
% run in folder with subfolders s1, s2, ..., s7
clear, close all
tic
folders = dir(); dirs = find(vertcat(folders.isdir));
folders = folders(dirs);

for foldi = 1:length(folders)
    direct = folders(foldi);
    if strcmp(direct.name(1), 's')
        cd (direct.name)
        sublist = dir('*_rsa.mat');
        sublist = {sublist.name};
        disp (['measurements -> ' num2str(length(sublist))]);
        for filei=1:length(sublist)
            load(sublist{filei});
            all{filei} = rsaZ;    
        end

        rsaZ = cat (2, all{:});

        cd .. % goes up one directory
        filename = [sublist{filei}(1:5) '_rsa'];
        save (filename, 'rsaZ');
    end
end


toc

%% LOAD r_H-LTC (enc or retrieval) PARTIAL COrrelation
clear, close all
tic
% Get a list of data files ready to be analyzed
sublist = dir('*_rsa.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
chan2plot = 1;

for subji=1:length(sublist)
    load(sublist{subji});
    %chan2plot(subji)
    %all{subji,1} = rsaZ ;
    all{subji,1} = squeeze(rsaZ(chan2plot, :, :, :)); 
    
end

timeBins1 = timeBins (:, [1, end]);

toc

%% VARIANCE ANALYSIS
% first create std_e_r_H std_e_r_LTC files from the all2all files 
%%load all2all
tic     
clearvars -except ALLEEG
% Get a list of data files ready to be analyzed
sublist = dir('*all2all.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subno=1:length(sublist)
    load(sublist{subno});
    s_all{subno} = all2all; 
    if exist ('combs')
        combs_all{subno} = combs;
    end
    if exist ('chanNames')
        chanNames_all{subno} = chanNames;
    end
end
toc

%% mean power per subject and condition all2all
clearvars -except s_all
ch = 2;
Nsubj = length(s_all);
for subji = 1:Nsubj
    all2all = s_all{subji};
    mE_perS(subji, :, :) = squeeze(mean(all2all(:,1,ch,:,:), 1));
    mR_perS(subji, :, :) = squeeze(mean(all2all(:,2,ch,:,:), 1));
    stdE_perS(subji, :, :) = squeeze(std(all2all(:,1,ch,:,:), [], 1));
    stdR_perS(subji, :, :) = squeeze(std(all2all(:,2,ch,:,:), [], 1));
    nTrials = size (all2all, 1); nPoints = size(all2all,5);
    allT_R = reshape (all2all, [2 2 44 nTrials*nPoints]);
    m_allT_R = squeeze(mean(allT_R(1, ch, :,:),4, 'omitnan'));
    std_allT_R = squeeze(std(allT_R(1, ch, :,:), [], 4, 'omitnan'));
    zS_R_perS(subji, :, :) = mean (bsxfun(@rdivide, (squeeze(all2all(:,2,ch,:,:)) ...
                        - mE_perS(subji, :, :)), stdR_perS(subji,:,:)));
%     zS_E_perS(subji, :, :) = ... 
%             mean (bsxfun(@rdivide, all2all(:,1,ch,:,:) - mean(all2all(:,1,ch,:,:),5, 'omitnan'), ...
%             std(all2all(:,1,ch,:,:),[], 5, 'omitnan')), 1);
%     zS_R_perS(subji, :, :) = ...
%             mean (bsxfun(@rdivide, all2all(:,2,ch,:,:) - mean(all2all(:,2,ch,:,:),5, 'omitnan'), ...
%             std(all2all(:,2,ch,:,:),[], 5, 'omitnan')), 1);        
end

mE = squeeze(mean(mE_perS));
mR = squeeze(mean(mR_perS));
stdE = squeeze(mean(stdE_perS));
stdR = squeeze(mean(stdR_perS));
%mZcE = squeeze(mean(zS_E_perS));
mZcR = squeeze(mean(zS_R_perS));


%% plot powerspectrum
min_freq1 = 1;max_freq1 = 29;num_frex1 = 29; min_freq2 = 30;max_freq2 = 100;num_frex2 = 15;baseline_window = [ -800 -300]; frex1 = linspace(min_freq1,max_freq1,num_frex1); frex2 = linspace(min_freq2,max_freq2,num_frex2); frex = [frex1 frex2]; num_frex = size(frex, 2); 
tftimes = (1:2:2251*2) - 500;

cond1 = mE;
cond2 = mR;

xlim = [-500 3000];
figure(1);
subplot (311)
contourf(tftimes,frex,cond1,40,'linecolor','none');hold on;
set(gca,'ydir','normal','xlim',xlim); %, 'clim',[-4 1]
%title (['Congruent - Encoding']);colorbar;
title (['Encoding']);colorbar;
subplot (312)
contourf(tftimes,frex,cond2,40,'linecolor','none');hold on;
set(gca,'ydir','normal','xlim',xlim); %'clim',[-15 15], 'clim',[-4 1]
%title (['Incongruent - Encoding']);colorbar;
title (['Retrieval']);colorbar;
subplot (313)
diff = cond1 - cond2;
contourf(tftimes,frex,diff,40,'linecolor','none');hold on;
set(gca,'ydir','normal','xlim',xlim, 'clim',[-3 3]); %'clim',[-15 15]
title (['Difference']);colorbar;

saveas (1,['SIDR_all.jpg']);





%% exclude subject 2
subj2exc = 2;
stdE_perS_H(subj2exc,:,:) = [];
stdR_perS_H(subj2exc,:,:) = [];
stdE_perS_LTC(subj2exc,:,:) = [];
stdR_perS_LTC(subj2exc,:,:) = [];


%% VARIANCE ANALYSIS
erTime      = 2000:4000;
b1          = 1:9;
b2          = 10:19;
b3          = 20:29;
b4          = 30:44;

var2useE = zS_R_perS_H; %
var2useR = zS_R_perS_LTC; %
% var2useE = stdR_perS_H; %
% var2useR = stdR_perS_LTC; %


%calculate variance of the specific band over time for each subject
std_b1_H    =   squeeze(mean (var2useE(:, b1,erTime), 2));
std_b2_H    =   squeeze(mean (var2useE(:, b2,erTime), 2));
std_b3_H    =   squeeze(mean (var2useE(:, b3,erTime), 2));
std_b4_H    =   squeeze(mean (var2useE(:, b4,erTime), 2));
m_std_b1_H    =   mean (std_b1_H,2);
m_std_b2_H    =   mean (std_b2_H,2);
m_std_b3_H    =   mean (std_b3_H,2);
m_std_b4_H    =   mean (std_b4_H,2);


std_b1_LTC    =   squeeze(mean (var2useR(:, b1,erTime), 2));
std_b2_LTC    =   squeeze(mean (var2useR(:, b2,erTime), 2));
std_b3_LTC    =   squeeze(mean (var2useR(:, b3,erTime), 2));
std_b4_LTC    =   squeeze(mean (var2useR(:, b4,erTime), 2));
m_std_b1_LTC    =   mean (std_b1_LTC,2);
m_std_b2_LTC    =   mean (std_b2_LTC,2);
m_std_b3_LTC    =   mean (std_b3_LTC,2);
m_std_b4_LTC    =   mean (std_b4_LTC,2);

%std_b1_LTC  =   stdR_perS_H(:, b1,erTime), 2);

%% plot  4 bar 
clear data
data.data = [m_std_b1_H m_std_b2_H m_std_b3_H m_std_b4_H];
%data.data = [m_std_b1_LTC m_std_b2_LTC m_std_b3_LTC m_std_b4_LTC];


figure(1); set(gcf,'Position', [0 0 400 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2 3 4], data.data, '.'); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',24);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'1-9 Hz  ',...
    '10-19 Hz  ', '20-29 Hz  ', '30-100Hz  '}, 'FontSize', ...
    18, 'linew',2, 'xlim',[0 5], 'ylim',[-.55 .35] );%[-0.17 0.1]
xticklabel_rotate;
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

export_fig(1, '_var_in_bands.png','-transparent', '-r300');
%close all;


%% plot  4 boxplots
clear data
data.data = [m_std_b1_H m_std_b2_H m_std_b3_H m_std_b4_H];
data.data = [m_std_b1_LTC m_std_b2_LTC m_std_b3_LTC m_std_b4_LTC];



figure(1); set(gcf,'Position', [0 0 400 400]); 
mean_S = mean(data.data, 1);
hb = boxplot (data.data); hold on;
%set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',24);hold on;
%set(hb,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'1-9 Hz  ',...
    '10-19 Hz  ', '20-29 Hz  ', '30-100Hz  '}, 'FontSize', ...
    18, 'linew',1, 'xlim',[0 5], 'ylim',[-0.17 0.1] );%[-2.5 2.5]
xticklabel_rotate;
plot(get(gca,'xlim'), [0 0],'k:','lineWidth', 1);

export_fig(1, '_var_in_bands.png','-transparent', '-r300');
%close all;
%% anova for variance

data4an = [dataH dataLTC];
%data4an = reshape(data4anova, [14, 4]);
data4anova = zeros(14,4);
data4anova(1:2:13,:)= dataH(1:1:7,:);
data4anova(2:2:14,:)= dataLTC(1:1:7,:);
[p,tbl] = anova2(data4anova,2);



%% compare 1-8 in H and LTC
close all; 
freq = 5; %2 = theta; 3 = alpha ..

H = data_H.data(freq,:); H = H';
LTC = data_LTC.data(freq,:); LTC = LTC';
data.data = [H LTC];

figure(1); set(gcf,'Position', [0 0 500 400]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',28);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LTC  '}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.01 0.1] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);
xticklabel_rotate;


[h p ci t] = ttest (data.data(:,1), data.data(:,2));

disp (['p = ' num2str(p) '   t= ' num2str(t.tstat)]);

filename = ['var' num2str(freq) '_.png'];
export_fig(1, filename,'-transparent', '-r300');
close all;

%%

data4anova  = data.data(:);
freqlab     = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
subjlab     = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7];
arealab     = [1:28];

% (Y,S,F1,F2,FACTNAMES) 
out_anova   = rm_anova2(data4anova, subjlab, freqlab);




%% calculate compare low/high ratio in H and LTC
erTime      = 1751:4000;
lfb         = 1:9;
hfb         = 10:44;
ratioRH = mean(stdR_perS_H(:, lfb,erTime), 2)   ./ mean(stdR_perS_H(:,hfb, erTime), 2);
ratioR_LTC = mean(stdR_perS_LTC(:, lfb,erTime), 2) ./ mean(stdR_perS_LTC(:,hfb, erTime), 2);
 
rH = squeeze(mean(ratioRH, 3));
rLTC = squeeze(mean(ratioR_LTC, 3));

[h p ci stats] = ttest(rH, rLTC);
p
stats.tstat


%% 2 bar
ylim = [.75 1.5];
%ylim = [-.1 .13]; %for the hc-lc comparison

data.data = [rH rLTC];
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(7);
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
hb = plot ([1 2 ], data.data(:,1:2)); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',25);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
%set(hb, 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LTC  '}, ...
    'FontSize', 25, 'linew',2, 'xlim', [0.4 2.6], 'ylim', ylim );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);
xticklabel_rotate;

set (gcf, 'Position', [200 200 400 450]);
export_fig(2, '_var_ratio.png','-transparent', '-r300');
close all;   


%% load contribution data / contribution of individual electrodes or frequencies
%this is done for each condition separately load electrodes data
clearvars 
tic

sublist = dir('*_rsa.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    load(sublist{subji});
    %chan2plot(subji)
    all{subji,1} = rsaZ_J_all; 
    %all{subji,1} = rsaZ_J;%when more than 1 electrode
    
end

cd .. % goes up one directory
filename = [sublist{subji}(5:end-13)];
%filename = [sublist{subji}(5:end-12)];
eval([filename '= all;']);
%save (filename, filename, '-v7.3');
save (filename, filename);

toc


%% determine most informative frequencies
clear mxE_all_minus mxE_all;
c_all = SI_all;
c_all_minus = SI_all_minus;

nSubj = length(c_all);
nChans = size(c_all_minus{1}, 1); %take first subject to determine 
nFreq = size(c_all_minus{1}, 2);  %nº of electrodes and frequencies
%mean across trials when including all and for each removal
for subji = 1:nSubj
    mxE_all_minus(subji, :, :) = mean(c_all_minus{subji}, 3);
    mxE_all(subji, :) = mean(c_all{subji}, 2);
end

clear p;
%contribution of each electrode for each subject
for subji = 1:nSubj
    for chani = 1:nChans
        for freqi = 1:nFreq
            p(subji, chani, freqi) = ...
                ( (mxE_all(subji,chani) - mxE_all_minus(subji, chani, freqi) ) ...
                                / mxE_all(subji,chani) ) * 100;
            %(mxE_all(subji,chani) - mxE_all_minus(subji, chani, freqi) );
        end
    end
end

p_H     = squeeze(p(:,1,:));
p_LTC   = squeeze(p(:,2,:));

%%zscore
p_LTCz = bsxfun(@rdivide, p_LTC - mean(p_LTC,2), std(p_LTC,[], 2));
p_Hz = bsxfun(@rdivide, p_H - mean(p_H,2), std(p_H,[], 2));


%%
figure();
subplot (211)
imagesc(p_Hz(:,1:44));hold on; colorbar;
subplot (212)
imagesc(p_LTCz(:,1:44));hold on; colorbar;

%% Frequency analysis in specific clusters
%load p_H p_LTC




%% remove subject
p_H(2,:) = []; p_LTC(2,:) = [];
p_Hz(2,:) = []; p_LTCz(2,:) = [];


%% 
d2plt_1 = p_Hz(:,1:44);
d2plt_2 = p_LTCz(:,1:44);

figure();
subplot (211);
plot(d2plt_1', 'Color', [.5 .5 .5]); hold on;
plot (mean(d2plt_1, 1), 'r', 'LineWidth', 2); hold on; 
%set (gca, 'ylim', [-2 5]);
subplot (212)
plot(d2plt_2', 'Color', [.5 .5 .5]); hold on;
plot (mean(d2plt_2, 1), 'r', 'LineWidth', 2);
%set (gca, 'ylim', [-2 5]);

%%
figure();
%diff = mean(d2plt_1) - mean(d2plt_2); 

mean1 = mean(d2plt_1, 1);
std_1 = std(d2plt_1, [], 1);
mean2 = mean(d2plt_2, 1);
std_2 = std(d2plt_2, [], 1);
se_1 = std_1 / sqrt(size (d2plt_1, 1));
se_2 = std_2 / sqrt(size (d2plt_2, 1));


%h1 = bar (mean1);hold on;
e1=errorbar(mean1, se_1,'r'); hold on;
%h2 = bar (mean2);hold on;
e2 = errorbar(mean2, se_2,'k'); hold on;

set (gca, 'xlim', [0.1 44]);
%hb = plot ([1 2 ], mean1); hold on; % > lines


[h p ci t] = ttest (d2plt_1, d2plt_2);

disp (['p = ' num2str(p) '    t = ' num2str(t.tstat)]);
p = p'; 

%%
plot (mean1, 'r', 'LineWidth', 1); hold on; 
plot (mean2, 'b', 'LineWidth', 1); hold on; 
%plot (diff, 'k', 'LineWidth', 1); hold on; 




%% rank electrodes according to their contribution
p(p==0) = NaN;
p_ranked = tiedrank(p); 
thres = prctile (p_ranked, 90); %gives the thres per column
idx90 = p_ranked>thres;

idxPositive = (p > 0);
numChan = sum (idxPositive)


%%
diffM = idxPositive_A - idxPositive_NA ;
am  = sum (diffM ==-1)
ana = sum (diffM ==1)
%diffM = idx90_A - idx90_NA ;
%%
am  = sum (idxPositive_NA)



%% find hipp cluster that correlates with LTC at encoding
test(test~=1.75)=0; test(test==1.75)=1;
test = flipud(test);
clustinfoHclust = bwconncomp(test); 
figure();
imagesc(test);axis square
numPixHclust = cellfun(@numel,clustinfoHclust.PixelIdxList);

%% test to understand the flipud thing

test = zeros (10);
test(2) = 1; test (12) = 1;
testfliped = flipud(test);
roiX = 1:3; roiY = 1:3;

mFp= test(roiX, roiY); mFp = mFp(:);mFp = mean (mFp);
mF = testfliped(roiX, roiY); mF = mF(:);mF = mean (mF);

%% exclude subjects from the saved contrasts
s = [2];
HC{s} = []; HC = HC(~cellfun('isempty',HC)); 
SISP{s} = []; SISP = SISP(~cellfun('isempty',SISP)); 
SIDR{s} = []; SIDR = SIDR(~cellfun('isempty',SIDR)); 

%% exclude subjects from the saved contrasts
load SI_H; load SI_LTC;
s = [2];
HC{s} = []; HC = HC(~cellfun('isempty',HC)); 
SI_H{s} = []; SI_H = SI_H(~cellfun('isempty',SI_H)); 
SI_LTC{s} = []; SI_LTC = SI_LTC(~cellfun('isempty',SI_LTC)); 
SISP_H{s} = []; SISP_H = SISP_H(~cellfun('isempty',SISP_H)); 
SISP_LTC{s} = []; SISP_LTC = SISP_LTC(~cellfun('isempty',SISP_LTC)); 



%% overlapping clusters 

clustInfoOverlapping = bwconncomp(test);



%% MEAN IN CLUSTER OF INTEREST
%load clustInfoReal_H ; load clustInfoReal_LTC_HC_LC;
clustInfoReal1 = clustInfoReal_LTC_SIDI;
clustInfoReal2 = clustInfoReal_LTC_HCLC;
%clustInfoReal = clustInfoReal_H;

rsa_cond1 = cell2mat(cellfun(@mean, HC, 'un',0)); 
rsa_cond2 = cell2mat(cellfun(@mean, LC, 'un',0)); 

takeclust = 1;
roiX  = 25:31; roiY  = 6:11; % % Hippocampus
%roiX = 25:31; roiY = 13:23; % LTC > real
pixel = 16;%H_cong_inc = 6; % LTC_HC_LC = 12 % LTC_SIDI = 16-19

%without baseline
% roiX = 20:26; roiY = 9:17; % LTC > 41-75 > use pixel= 8
% pixel = 8;

mlim = 36:75;
rsa_cond1 = rsa_cond1(:,mlim,mlim);
rsa_cond2 = rsa_cond2(:,mlim,mlim);

nSubj = size(rsa_cond1, 1);
if takeclust
    ROI_cond1 = rsa_cond1(:,clustInfoReal.PixelIdxList{pixel}); 
    ROI_cond2 = rsa_cond2(:,clustInfoReal.PixelIdxList{pixel}); 
    m_roi_cond1 = mean(ROI_cond1, 2);
    m_roi_cond2 = mean(ROI_cond2, 2);
else
    ROI_cond1 = rsa_cond1(:,roiX, roiY); 
    ROI_cond2 = rsa_cond2(:,roiX, roiY); 
    ROI_cond1 = reshape (ROI_cond1, [nSubj size(roiX,2) * size(roiY, 2)]);
    ROI_cond2 = reshape (ROI_cond2, [nSubj size(roiX,2) * size(roiY, 2)]);
    m_roi_cond1 = mean(ROI_cond1, 2);
    m_roi_cond2 = mean(ROI_cond2, 2);
end

test = squeeze(mean(rsa_cond1, 1));
test(1:end) = 0;
%test(roiX,roiY,:) = .5;
test(clustInfoReal1.PixelIdxList{16}) = 1;
test(clustInfoReal2.PixelIdxList{12}) = test(clustInfoReal2.PixelIdxList{12}) +1;
test(test ~=2) = 0;
test = flipud(test);
ids = find (test == 2);


figure(1);imagesc(test);axis square;hold on;colorbar;
%figure(1);contourf(test, 40, 'linecolor', 'none');axis square;hold on;
% plot([mfL+0.5 mfL+0.5],get(gca,'ylim'),'w', 'LineWidth', 2); 
% plot(get(gca,'xlim'), [bins - (mfL-0.5) bins - (mfL-0.5)],'w', 'LineWidth', 2);
% plot(get(gca,'xlim'), [bins+.5 .5],'w', 'LineWidth', 2); % diagonal
% set (gca, 'clim', [0 1]);
%remove labels
axesHandles = findall(0, 'type', 'axes');
for i=1:length(axesHandles)
   set (axesHandles(i), 'visible', 'off'); % remove labels 
end


[h p ci t] = ttest (m_roi_cond2);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set (gcf, 'Position', [200 200 400 450]);
export_fig(2, 'plotting1.png','-transparent', '-r300');
%close all;   

%% PLOT 2 bar
%ylim = [-.05 .1];
ylim = [-.125 .15]; %for the hc-lc comparison

data.data = [m_roi_cond1 m_roi_cond2];
figure(2); set(gcf,'Position', [0 0 560 500]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(7);
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
hb = plot ([1 2 ], data.data(:,1:2)); hold on; % > lines
set(hb, 'lineWidth', 2, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
%set(hb,'linestyle','none', 'lineWidth', 2);
%set(hb, 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'HC  ', 'LC  '}, ...
    'FontSize', 32, 'linew',3, 'xlim', [0.4 2.6], 'ylim', ylim );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);
xticklabel_rotate;

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

export_fig(2, '_currentBarPlot_m2C.png','-transparent', '-r300');
close all;   




%% predictive analysis
clearvars
load carsmall

tbl = table(Weight,Acceleration,MPG,'VariableNames',{'Weight','Acceleration','MPG'});



lm = fitlm(tbl,'MPG~Weight+Acceleration')





%% %%%%> RECOGNITION  
%clear, close all
clearvars -except ALLEEG 

chan2plot = [1 1 1 1 1 1 1 1];
%chan2plot = [2 2 2 2 2 2 2 2 ];



tic
% Get a list of data files ready to be analyzed
sublist = dir('*_rsa.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    load(sublist{subji});
    %chan2plot(subji)
    all{subji,1} = squeeze(rsaZ (chan2plot(subji),:,:,:)); %when more than 1 electrode
    %all{subji,1} = rsaZ;
    
end

if exist ('timeBins')
    timeBins1 = timeBins (:, [1, end]);
end

cd .. % goes up one directory
if strcmp (sublist{subji}(6), '_')
    disp ('short name');
    %filename = [sublist{subji}(4:end-14)]; % es 12 o 14
    %filename = [sublist{subji}(4:end-36)]; % es 12 o 14
    %filename = [sublist{subji}(4:end-8)]; % for DI batch files
    filename = [sublist{subji}(4:end-34)]; % 
else
    disp ('long name');
    filename = [sublist{subji}(4:end-34)]; % es 12 o 14
    %filename = [sublist{subji}(4:end-12)]; % es 12 o 14
end
eval([filename '= all;']);
save (filename, filename, '-v7.3');
%save (filename, filename);


toc








%%