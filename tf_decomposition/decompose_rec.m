
function [oneListPow] = decompose_rec (oneListTraces_c, cfg_dec)
        
timeRes     =       cfg_dec.timeRes; 
blne        =       cfg_dec.blne;    
finalCut    =       cfg_dec.finalCut;
   
sr = 500;
data_ft = mat2ft(oneListTraces_c, sr);
        
size(oneListTraces_c)
if strcmp (timeRes, 'all')
    disp ('Final Cut: '); 
    ble1    = floor ((size(oneListTraces_c, 2) / 2) +  blne(1) * 500)
    ble2    = floor ((size(oneListTraces_c, 2) / 2) +  blne(2) * 500)
    lim_1   = floor ((size(oneListTraces_c, 2) / 2) +  finalCut(1) * 500)
    lim_2   = floor ((size(oneListTraces_c, 2) / 2) +  finalCut(2) * 500)
end

if timeRes == 0.01
    disp ('10ms'); 
    total = (size(oneListTraces_c, 2) * 2 ) / 10; %-6-6 > 1200
    ble1    = floor (total/2 +  blne(1)*100 )
    ble2    = floor (total/2 +  blne(2)*100)
    lim_1   = floor (total/2 +  finalCut(1)*100) +1;
    lim_2   = floor (total/2 +  finalCut(2)*100)
end

if timeRes == 0.05
    disp ('50ms');  
    total = (size(oneListTraces_c, 2) * 2 ) / 50 %-6-6 > 240
    ble1    = floor (total/2 +  blne(1)*20 )
    ble2    = floor (total/2 +  blne(2)*20)
    lim_1   = floor (total/2 +  finalCut(1)*20)
    lim_2   = floor (total/2 +  finalCut(2)*20)
end


disp ('>> extracting power ... ');

%fieldtrip analysis
cfg              = [];
cfg.method       = 'wavelet'; %%'mtmconvol'; % 'wavelet'; %
cfg.width        = linspace(3, 6, 29);
cfg.output       = 'pow';
cfg.foi          = [1:1:29];   % analysis 2 to X Hz in steps of 2 Hz 
if strcmp (timeRes, 'all')
    cfg.toi          = 'all'; % takes as reference the number of time windows defined above
else
    cfg.toi          = 0:timeRes:(length(oneListTraces_c)/sr); 
end
cfg.keeptrials   = 'yes'; % keep individual trials or average 
tf_data_L          = ft_freqanalysis(cfg, data_ft);
dataL = tf_data_L.powspctrm;

size(dataL)


cfg              = [];
cfg.method       = 'wavelet'; %'mtmconvol'; % 'wavelet'; %
cfg.width        = linspace(6, 12, 15); 
cfg.output       = 'pow';
cfg.foi          = [30:5:100];   
if strcmp (timeRes, 'all')
    cfg.toi          = 'all'; 
else
    cfg.toi          = 0:timeRes:(length(oneListTraces_c)/sr); % 50ms; 
end
cfg.keeptrials   = 'yes'; % keep individual trials or average 
tf_data_H          = ft_freqanalysis(cfg, data_ft);
dataH = tf_data_H.powspctrm;


size(dataH)


% decibel normalization
dataL =  10*log10( bsxfun(@rdivide, dataL, mean(dataL(:,:,:,ble1:ble2),4, 'omitnan') )); %decibel norm   
dataH =  10*log10( bsxfun(@rdivide, dataH, mean(dataH(:,:,:,ble1:ble2),4, 'omitnan') )); %decibel norm   

dataLH = cat (3, dataL, dataH);


%cut edge artefacts from -6-6 to -3-3
dataLH = dataLH(:,:,:,lim_1:lim_2);
oneListPow = dataLH;















        

