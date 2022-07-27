
function [] = calculate_rsa(files, win_width, mf, onebyone)

if onebyone
    s_all           = files.s_all; 
    chanNames_all   = files.chanNames_all;
    sublist         = files.sublist;

    %%Vector version (used in Recognition)
    datetime
    win_width           =  win_width; 
    mf                  =  mf; 
    f                   =  [1:44];
    surrogates          =  0; %generates surrogates for the DI conditoin (activate only with SI)
    initial_n2s         =  0; %if zero takes max possible for each max 500 because of S7
    fact4dwsp           =  1; %downsample factor for feature vector (1 = no downsampling)

    disp (['>>> win of ' num2str(win_width*2) 'ms sliding in ' num2str(mf*2) 'ms steps... ']);
    tic;
    nTimepoints         = size (s_all{1}, 5);
    bins                =  floor ( (nTimepoints/mf)- win_width/mf+1 );
    timeBins            = zeros (bins,win_width);
    Nsubj               = length(s_all);
    nFreq               = length(f);

    for subji = 1:Nsubj 
        fprintf('\n');
        all2all = s_all{subji};
        n2s     = initial_n2s; %reset number of surrogates
        if ndims(all2all) == 5 %excludes all2all files with no trials
            disp (['File ' num2str(subji) '  ' num2str(size(all2all, 3)) ' electrodes  ' ...
                , num2str(size(all2all, 1)) ' trials']);

            if ~surrogates 
                disp (['Subj ' num2str(subji) '/' num2str(Nsubj) '    ' num2str(size(all2all, 3)) ' electrodes  '])
                chanNames = chanNames_all{subji};
                trialN = size(all2all, 1);
                chanN = size(all2all, 3);
                x_prev = zeros (1, 1, 1, nFreq, win_width); 
                y_prev = zeros (1, 1, 1, nFreq, win_width); 
                x = zeros (1, win_width * nFreq);
                y = zeros (1, win_width * nFreq);
                xM = zeros (chanN, trialN,(win_width* nFreq) / fact4dwsp, bins);
                yM = zeros (chanN, trialN, (win_width* nFreq) / fact4dwsp, bins);
                mX = zeros (win_width * nFreq, bins);
                mY = zeros (win_width * nFreq, bins);
            else
                disp ('Surrogates');
                disp (['Subj ' num2str(subji) '/' num2str(Nsubj) '    ' num2str(size(all2all, 3)) ' electrodes  '])
                %t1          = 1:size(all2all, 1); % no repetitions
                %combs       = nchoosek(t1, 2);    % no repetitions
                t1          = size(all2all, 1); %repetitions are needed because 
                combs       = combinator(t1, 2, 'p'); %each item has enc and retrival
                if n2s == 0 n2s = length(combs); end
                if n2s >= length(combs) n2s = length(combs); else n2s = initial_n2s; end
                selCombs    = randsample(length(combs), n2s);
                chanN = size(all2all, 3);
                x_prev = zeros (1, 1, 1, nFreq, win_width); 
                y_prev = zeros (1, 1, 1, nFreq, win_width); 
                x = zeros (1, win_width * nFreq);
                y = zeros (1, win_width * nFreq);
                trialN = length(selCombs);
                xM = zeros (chanN, trialN, (win_width* nFreq) / fact4dwsp, bins);
                yM = zeros (chanN, trialN,  (win_width* nFreq) / fact4dwsp, bins);
                mX = zeros (win_width * nFreq, bins);
                mY = zeros (win_width * nFreq, bins);

            end


            if ~surrogates
                for chani = 1:size(all2all, 3) 
                    fprintf('\n');fprintf('Channel:      ');
                    fprintf('%d %s', chani, ' '); if (rem(chani,20) == 0) fprintf('\n'); end
                    fprintf('\n'); fprintf('Trial epoching:      '); %check this later
                for triali = 1:trialN %trials
                    if (triali < 10) fprintf('\b\b\b'); elseif (triali < 100) fprintf('\b\b\b\b'); else fprintf('\b\b\b\b\b');  end
                        fprintf('%d %s', triali, ' '); %if (rem(triali,20) == 0) fprintf('\n'); end
                        for timei = 1:bins 
                            timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                            x_prev = all2all(triali, 1,chani,f,timeBins(timei,:));
                            x = x_prev(:);
                            %x = downsample(x, fact4dwsp); 
                            xM(chani, triali, :, timei) =  x;
                            y_prev = all2all(triali, 2,chani,f,timeBins(timei,:));
                            y = y_prev(:);
                            %y = downsample (y, fact4dwsp); 
                            yM(chani, triali, :, timei) =  y; % so I don0t have to transpose later
                        end
                    end
                end
                fprintf('\n');

                rsaZ = zeros (chanN, trialN, bins, bins);
                for chi = 1:chanN
                    fprintf('\n'); fprintf('Channel:      '); %check this later
                    fprintf('%d %s', chi, ' '); if (rem(chi,20) == 0) fprintf('\n'); end
                    fprintf('\n'); fprintf('trial correlation:      '); %check this later
                    for ti = 1:trialN
                        if (ti < 10) fprintf('\b\b\b'); elseif (ti < 100) fprintf('\b\b\b\b'); 
                        else fprintf('\b\b\b\b\b'); end
                        fprintf('%d %s', ti, ' '); %if (rem(triali,20) == 0) fprintf('\n'); end
                        mX= squeeze(xM(chi, ti,:,:));
                        mY= squeeze(yM(chi, ti,:,:));
                        [r p] = corr (mX, mY,'Type','s'); %Spearman
                        rsaZ(chi, ti, :, :) = atanh(r);
                    end
                end

            else %below for surrogates

                for chani = 1:size(all2all, 3) 
                    fprintf('\n');fprintf('Channel:      ');
                    fprintf('%d %s', chani, ' '); if (rem(chani,20) == 0) fprintf('\n'); end
                    fprintf('\n'); fprintf('Trial epoching:      '); %check this later
                    for triali = 1:length(selCombs) %trials
                        if (triali < 10) fprintf('\b\b\b'); elseif (triali < 100) fprintf('\b\b\b\b'); else fprintf('\b\b\b\b\b');  end
                            fprintf('%d %s', triali, ' '); %if (rem(triali,20) == 0) fprintf('\n'); end
                            for timei = 1:bins 
                                timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                                x_prev = all2all(combs(selCombs(triali), 1), 1,chani,f,timeBins(timei,:));
                                x = x_prev(:);
                                %x = downsample (x, fact4dwsp); 
                                xM(chani, triali, :, timei) =  x;
                                y_prev = all2all(combs(selCombs(triali), 2), 2,chani,f,timeBins(timei,:));
                                y = y_prev(:);
                                %y = downsample (y, fact4dwsp); 
                                yM(chani, triali, :, timei) =  y; % so I don0t have to transpose later
                            end
                        end
                    end
                    fprintf('\n');

                    rsaZ = zeros (chanN, trialN, bins, bins);
                    for chi = 1:chanN
                        fprintf('\n'); fprintf('Channel:      '); %check this later
                        fprintf('%d %s', chi, ' '); if (rem(chi,20) == 0) fprintf('\n'); end
                        fprintf('\n'); fprintf('trial correlation:      '); %check this later
                        for ti = 1:trialN
                            if (ti < 10) fprintf('\b\b\b'); elseif (ti < 100) fprintf('\b\b\b\b'); 
                            else fprintf('\b\b\b\b\b'); end
                            fprintf('%d %s', ti, ' '); %if (rem(triali,20) == 0) fprintf('\n'); end
                            mX= squeeze(xM(chi, ti,:,:));
                            mY= squeeze(yM(chi, ti,:,:));
                            [r p] = corr (mX, mY,'Type','s'); %Spearman
                            rsaZ(chi, ti, :, :) = atanh(r);
                        end
                    end
            end

            if ~surrogates
                fname = strsplit(sublist{subji}(1:end), '_');
                if strcmp(fname{2}, 'DI') || strcmp(fname{2}, 'DIC') 
                    filename = [fname{1} '_' fname{2} '_' fname{3}   '_rsa.mat'];
                elseif strcmp(fname{2}, 'DIN') || strcmp(fname{2}, 'DINC') || ...
                    strcmp(fname{2}, 'DINHCC') || strcmp(fname{2}, 'SIDRDINC') 
                    filename = [fname{1} '_' fname{2} '_' fname{4}   '_rsa.mat'];
                else
                filename = [fname{1} '_' fname{2} '_rsa.mat'];
            end
            else
                fname = strsplit(sublist{subji}(1:end), '_');
                filename = [fname{1} '_' strrep(fname{2}, 'S', 'D') '_OBO'   '_rsa.mat'];
            end 
            save (filename, 'rsaZ', 'chanNames', 'timeBins'); 
            fprintf('\n');
        end
    end
    toc
    
    
else %global RSA for channels included 
    
    disp ('global RSA');
    s_all           = files.s_all; 
    chanNames_all   = files.chanNames_all;
    sublist         = files.sublist;

    %%Vector version (used in Recognition)
    win_width           =  win_width; 
    mf                  =  mf; 
    f                   =  [1:44];
    surrogates          =  0; %generates surrogates for the DI conditoin (activate only with SI)
    initial_n2s         =  0; %if zero takes max possible for each max 500 because of S7
    fact4dwsp           =  1; %downsample factor for feature vector (1 = no downsampling)

    disp (['>>> win of ' num2str(win_width*2) 'ms sliding in ' num2str(mf*2) 'ms steps... ']);
    tic;
    nTimepoints         = size (s_all{1}, 5);
    bins                =  floor ( (nTimepoints/mf)- win_width/mf+1 );
    timeBins            = zeros (bins,win_width);
    Nsubj               = length(s_all);
    nFreq               = length(f);

    for subji = 1:Nsubj 
        fprintf('\n');
        all2all = s_all{subji};
        n2s     = initial_n2s; %reset number of surrogates
        disp (['File ' num2str(subji) '  ' num2str(size(all2all, 3)) ' electrodes  ' ...
            , num2str(size(all2all, 1)) ' trials']);
        disp (['Subj ' num2str(subji) '/' num2str(Nsubj) '    ' num2str(size(all2all, 3)) ' electrodes  '])
        chanNames = chanNames_all{subji};
        trialN = size(all2all, 1);
        nChan = size(all2all, 3);
        x_prev = zeros (1, 1, 1, nFreq, win_width); 
        y_prev = zeros (1, 1, 1, nFreq, win_width); 
        x = zeros (1, win_width * nFreq);
        y = zeros (1, win_width * nFreq);
        xM = zeros (trialN,(win_width* nFreq *nChan) / fact4dwsp, bins);
        yM = zeros (trialN, (win_width* nFreq*nChan) / fact4dwsp, bins);
        mX = zeros (win_width * nFreq, bins);
        mY = zeros (win_width * nFreq, bins);

        fprintf('\n'); fprintf('Epoching:      '); %check this later
        for triali = 1:trialN %trials
            if (triali < 10) fprintf('\b\b\b'); elseif (triali < 100) fprintf('\b\b\b\b'); else fprintf('\b\b\b\b\b');  end
                fprintf('%d %s', triali, ' '); %if (rem(triali,20) == 0) fprintf('\n'); end
            for timei = 1:bins 
                timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                x_prev = all2all(triali, 1,:,f,timeBins(timei,:));
                x = x_prev(:);
                xM(triali, :, timei) =  x;
                y_prev = all2all(triali, 2,:,f,timeBins(timei,:));
                y = y_prev(:);
                yM(triali, :, timei) =  y; % so I don0t have to transpose later
            end
        end
        
        rsaZ = zeros (trialN, bins, bins);
        fprintf('\n');
        fprintf('\n'); fprintf('Correlation:      '); %check this later
        for ti = 1:trialN
            if (ti < 10) fprintf('\b\b\b'); elseif (ti < 100) fprintf('\b\b\b\b'); 
            else fprintf('\b\b\b\b\b'); end
            fprintf('%d %s', ti, ' '); %if (rem(triali,20) == 0) fprintf('\n'); end
            mX= squeeze(xM(ti,:,:));
            mY= squeeze(yM(ti,:,:));
            [r p] = corr (mX, mY,'Type','s'); %Spearman
            rsaZ(ti, :, :) = atanh(r);
        end



        fname = strsplit(sublist{subji}(1:end), '_');
        if strcmp(fname{2}, 'DI') || strcmp(fname{2}, 'DIC') 
            filename = [fname{1} '_' fname{2} '_' fname{3}   '_rsa.mat'];
        elseif strcmp(fname{2}, 'DIN') || strcmp(fname{2}, 'DINC') || strcmp(fname{2}, 'DINHCC') 
            filename = [fname{1} '_' fname{2} '_' fname{4}   '_rsa.mat'];
        else
            filename = [fname{1} '_' fname{2} '_rsa.mat'];
        end

        save (filename, 'rsaZ', 'chanNames', 'timeBins'); 
        fprintf('\n');
        end
    end
    toc
   
    
    
    
    
end
