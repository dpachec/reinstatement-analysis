
function [ALLEEG allR] = load_all_data_R (cfg)

only1sub        =   cfg.only1subj;
align2resp      =   cfg.align2resp; 
nIQR            =   cfg.nIQR; 
plotTrajEvents  =   cfg.plotTrajEvents;

%%original
ch2p{1}     = {'B1'; 'B8'}; 
ch2p{2}     = {'B''1'; 'B''7'}; 
ch2p{3}     = {'B''1'; 'B''7'}; 
ch2p{4}     = {'B''1'; 'B''9'}; 
ch2p{5}     = {'B''1'; 'B''9'}; 
ch2p{6}     = {'B1'; 'B8'}; 
ch2p{7}     = {'B1'; 'B9'}; 
ch2p{8}     = {'B''1'; 'B''10'}; 
ch2p{9}     = {'B1'; 'B10'}; 
ch2p{10}    = {'B''1'; 'B''6'}; 
ch2p{11}    = {'B''1'; 'B''7'}; 

% trials 2 exc were lost and recovered (computer stolen until trial 80)
% results involving novel items (tr > 80) are not included (except in s11)
tr2exc{1}   = [2 8 15 17 21 25 ]; % PAPER VERSION
tr2exc{2}   = [];
tr2exc{3}   = [1 2 3 18 19 29 35 39 ]; %
tr2exc{4}   = [6 7 8 10 12 13 14 15 16 19 20 21 26 27 32 36 37 39 40]; 
tr2exc{5}   = [];
tr2exc{6}   = [2 9 20 32 33 35 37 40]; % PAPER 
tr2exc{7}   = [];
tr2exc{8}   = [5 9 11 12 14 15 21 22 23 25 27 28 34 36 39]; 
tr2exc{9}   = []; %no exclusions from visual inspection
tr2exc{10}  = []; %
tr2exc{11}  = [1 2 3 6 9 10 16 17 19 20 21 29 31 33 34 35 36 86 117];


% tr2exc{1}   = [2 8 15 17 21 25 41 43 50 53 64 65]; % 
% tr2exc{2}   = [];
% tr2exc{3}   = [1 2 3 4 9 13 18 19 29 35 39 48 50 54 56 60 68 78]; %
% tr2exc{4}   = [3 4 6 7 8 10 12 13 14 15 16 17 19 20 21 26 27 32 33 36 37 39 40 48 51 52 53 ...
%                  54 55 56 57 60 63 69 72 74 76 78 80]; %
% tr2exc{5}   = [];
% tr2exc{6}   = [2 5 6 9 16 20 25 32 33 35 37 40 41 43 54 62 65 74 76 80]; % 
% tr2exc{7}   = [];
% tr2exc{8}   = [5 7 12 14 15 21 22 23 27 28 36 39 50 58 ...
%                63 65 70 74 85 95 97 99 100]; 
% tr2exc{9}   = []; 
% tr2exc{10}  = []; 
% tr2exc{11}  = [1 18 19 20 23 25 26 29 34 39 45 48 57 58 61 68 83 86 89 90 92 96 105 107 ]; 





%%Load S1 
subj = 's01';
findCvalue = 1;
offset = 0; 
tEr = 739.2;
strEvent = 'X < -900';
datafile_name = 'ASJ.set';
datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s1_log.txt';
datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s1_traj.txt';
chans2exc = [{'Evento';'TTL';'C127';'C128'; 'T4'; 'T5'; 'D3'; 'D4'; 'D5'; 'D6'; 'E3'; 'E4'; 'E5'; 'E6'; 'J15'}];
chans2plot = ch2p{1}; 
epo2exc = tr2exc{1};
%epo2exc = []; 

[EEG1 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                            tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
EEG1.subject        = subj; 
EEG1.epo2exc        = epo2exc; 
EEG1.chans2exc      = chans2exc;
EEG1.chans2plot     = chans2plot;
EEG1.allR           = allR;

if ~only1sub
    %%Load S2 
    subj = 's02'; 
    findCvalue = 1;
    offset = 0; 
    tEr = 1;
    strEvent = 'X > 500';
    datafile_name = 'MMM.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s2_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s2_traj.txt';
    chans2exc = [{'TTL'; 'Evento';'ECG D';'ECG I';'C220';'C221'; 'P''1';'P''2';'P''3';'P''4'; 'O''1';'O''2';'O''3';'O''4';'O''5';'O''6';'O''7'; 'J''1';'J''2';'J''3'; 'K''1';'K''2';'K''3'}];
    %chans2exc = [{'TTL'; 'Evento';'ECG D';'ECG I';'C220';'C221';'P''1';'P''2';'P''3';'P''4'; 'O''1';'O''2';'O''3';'O''4';'O''5';'O''6';'O''7'; 'J''1';'J''2';'J''3'; 'K''1';'K''2';'K''3';'R''1';'R''2';'R''3';'R''4';'R''5';'R''6';'R''7';'R''8';'R''9';'R''10';'R''11';'R''12';'R''14';'R''15';'S''1';'S''4';'S''6';'S''7';'S''8';'P''1';'P''3';'P''4';'P''6';'P''7';'M''2';'M''3';'M''7';'M''8';'M''9';'M''13';'M''15';'J''1';'J''6';'J''8';'J''9';'J''12';'K''1';'K''6';'K''7';'K''8';'K''9';'K''11';'K''12';'L''2';'L''3';'L''4';'L''5';'L''8';'B''2';'B''3';'B''4';'B''5';'B''6';'B''7';'B''8';'B''10';'C''1';'C''4';'C''7';'O''4';'O''7';'O''8';'O''13';'O''15';'Q''2';'Q''4';'Q''7';'Q''8';'Q''11';'S2';'S3';'S4';'S5';'S7';'S8';'S10';'S11';'S12';'P2';'P9';'P10';'K1';'K7';'K8';'L1';'L5';'L6';'L7';'L8';'J2';'J3';'J8';'J9';'J13';'J14';'J15';'C7';'Q3';'Q6';'Q8';'Q13';'O3';'O4';'O5';'O6';'O7';'O8';'O9'}];
    chans2plot = ch2p{2};% {'B''1'; 'B''7'}; 
    epo2exc = tr2exc{2};

    [EEG2 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                          tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG2.subject        = subj; 
    EEG2.epo2exc        = epo2exc; 
    EEG2.chans2exc      = chans2exc;
    EEG2.chans2plot     = chans2plot;
    EEG2.allR           = allR;

    %%Load S3
    subj = 's03'; 
    findCvalue = 1;
    offset = 93.5*500;
    tEr = 4807.3;
    strEvent = 'X < - 700';
    datafile_name = 'RGE.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s3_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s3_traj.txt';
    chans2exc = [{'Evento';'TTL';'C127';'C128'; 'T''10';'T''11';'T''12';'T''13';'T''14'; 'E''6';'E''7';'E''8'; 'E''9'}];
    chans2plot = ch2p{3};
    epo2exc = tr2exc{3};
   

    [EEG3 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                          tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG3.subject        = subj; 
    EEG3.epo2exc        = epo2exc; 
    EEG3.chans2exc      = chans2exc;
    EEG3.chans2plot     = chans2plot;
    EEG3.allR           = allR;
    EEG3.H              = 'C''1';

    %%Load S4
    subj = 's04'; 
    findCvalue = 1;
    offset = 0;
    tEr = 146.3;
    strEvent = 'X > 700';
    datafile_name = 'RCM.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s4_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s4_traj.txt';
    chans2exc = [{'Evento';'C128'; 'TTL';'C138';'C139';'C140';'C141';'C142';'ECG1';'ECG2';'C145'; 'A''1';'A''2';'A''3'}];
    chans2plot = ch2p{4};
    epo2exc = tr2exc{4};
    [EEG4 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                          tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG4.subject        = subj; 
    EEG4.epo2exc        = epo2exc; 
    EEG4.chans2exc      = chans2exc;
    EEG4.chans2plot     = chans2plot;
    EEG4.allR           = allR;

    %%Load S5
    subj = 's05'; 
    findCvalue = 1;
    offset = 0;
    tEr = 2397.5;
    strEvent = 'X > 900';
    datafile_name = 'BRM.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s5_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s5_traj.txt';
    chans2exc = [{'Evento';'TTL';'C123';'C124';'C125';'C126';'C127';'C128'; 'L''7';'L''8';'L''9';'L''10';'L''11'; 'L''12'; 'K''1';'K''2';'K''3';'K''4';'K''5';'K''6';'K''7';'K''8';'K''9';'K''10';'K''11';'K''12';'K''13';'K''14'; 'K''15'}];
    chans2plot = ch2p{5}; 
    epo2exc = tr2exc{5}; 
    
    [EEG5 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                          tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG5.subject        = subj; 
    EEG5.epo2exc        = epo2exc; 
    EEG5.chans2exc      = chans2exc;
    EEG5.chans2plot     = chans2plot;
    EEG5.allR           = allR;

    %%Load S6
    subj = 's06'; 
    findCvalue = 1;
    offset = 0;
    tEr = 72;
    strEvent = 'X > 900';
    datafile_name = 'BBB.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s6_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s6_traj.txt';
    chans2exc = [{'TTL'; 'V''9'; 'Evento';'C169';'C170';'C171'; 'ECG1';'ECG2';'C64'; 'T''1';'T''2';'T''3';'T''4';'T''5';'T''6';'T''7';'T''8';'T''9';'T''10';'T''11';'T''12'; 'A''7';'A''8';'A''9';'A''10';'A''11';'E''4';'E''5';'E''6';'E''7';'E''8'}];
    chans2plot = ch2p{6};
    epo2exc = tr2exc{6};
    
    [EEG6 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                          tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG6.subject        = subj; 
    EEG6.epo2exc        = epo2exc;
    EEG6.chans2exc      = chans2exc;
    EEG6.chans2plot     = chans2plot;
    EEG6.allR           = allR;


    %%Load S7
    subj = 's07'; 
    findCvalue = 1;
    offset = 0;
    tEr = 690;
    strEvent = 'X > 500';
    datafile_name = 'CMM.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s7_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s7_traj.txt';
    chans2exc = {'ECG1';'ECG2';'TTL';'C125';'C126';'C127';'C128'};
    chans2plot = ch2p{7};
    epo2exc = tr2exc{7};

    [EEG7 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                            tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG7.subject        = subj; 
    EEG7.epo2exc        = epo2exc;
    EEG7.chans2exc      = chans2exc;
    EEG7.chans2plot     = chans2plot;
    EEG7.allR           = allR;

    
    %%Load S8
    subj = 's08'; 
    offset = 1; 
    findCvalue = 1;
    tEr = 732; %multiplied by sampling rate in function
    strEvent = 'X < -900';
    datafile_name = 'CQJ.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s8_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s8_traj.txt';
    chans2exc = {'Evento'; 'Fz';'Cz';'Pz';'Fp1';'F7';'T3';'T1';'T5';'O1';'F3';'C3';'P3';'Fp2';'F8';'T4';'T2';'T6';'O2';'F4';'C4';'P4';'C126';'ECG1';'ECG2'};
    chans2plot = ch2p{8};
    epo2exc = tr2exc{8};
    
    [EEG8 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                                tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG8.subject        = subj; 
    EEG8.epo2exc        = epo2exc; 
    EEG8.chans2exc      = chans2exc;
    EEG8.chans2plot     = chans2plot;
    EEG8.allR           = allR;
    
    %%Load S9
    subj = 's09'; 
    findCvalue = 1;
    offset = -600; 
    tEr = 1; %multiplied by sampling rate in function
    strEvent = 'X > 500';
    datafile_name = 'VML.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s9_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s9_traj.txt';
    chans2exc = [{'Evento';'TTL'}];
    chans2plot = ch2p{9};
    epo2exc = tr2exc{9};
    [EEG9 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                                tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG9.subject        = subj; 
    EEG9.epo2exc        = epo2exc; 
    EEG9.chans2exc      = chans2exc;
    EEG9.chans2plot     = chans2plot;
    EEG9.allR           = allR;
    
    %%Load S10
    subj = 's10'; 
    findCvalue = 1;
    offset = 18.3 * 500; 
    tEr = 1; %multiplied by sampling rate in function
    strEvent = 'X < -800';
    datafile_name = 'CGA.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s10_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s10_traj.txt';
    chans2exc = [{'ECG1';'ECG2'; 'Fpz';'Fp1';'AFz';'AF3';'AF7';'Fz';'F1';'F3';'F5';'F7';'FCz';'FC1';'FC3';'FC5';'FT7';'Cz';'C3';'C5';'T3';'C4';'ROC';'LOC';'Chin1';'Chin2';'TTL'}];
    chans2plot = ch2p{10};
    epo2exc = tr2exc{10};
    [EEG10 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                                tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG10.subject        = subj; 
    EEG10.epo2exc        = epo2exc; 
    EEG10.chans2exc      = chans2exc;
    EEG10.chans2plot     = chans2plot;
    EEG10.allR           = allR;
    
    %%Load S11
    subj = 's11'; 
    findCvalue = 1;
    offset = 1; 
    tEr = 236.2; %multiplied by sampling rate in function
    strEvent = 'X > 500';
    datafile_name = 'ZZM.set';
    datafile_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionLFP';
    datalog_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s11_log.txt';
    datatraj_path = 'C:\Users\Neuropsychology\Documents\HOSPITAL\data\recognition\recognitionBehavioural\patients\data_all\s11_traj.txt';
    chans2exc = [{'DC6'}];
    chans2plot = ch2p{11};
    epo2exc = tr2exc{11};
    [EEG11 allR]         = ld_R (datafile_name, datafile_path, datalog_path, datatraj_path, ...
                                tEr, offset, align2resp, strEvent, subj, plotTrajEvents, findCvalue);
    EEG11.subject        = subj; 
    EEG11.epo2exc        = epo2exc; 
    EEG11.chans2exc      = chans2exc;
    EEG11.chans2plot     = chans2plot;
    EEG11.allR           = allR;

    ALLEEG = {EEG1 EEG2 EEG3 EEG4 EEG5 EEG6 EEG7 EEG8 EEG9 EEG10 EEG11}';

else
    ALLEEG = {EEG1}';
end

disp ([newline 'raw data loaded' string(datetime) newline]);


