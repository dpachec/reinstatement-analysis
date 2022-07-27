function [rsaZ] =  rsa_m(all2all, win_width, mf, f, meanInTime)


nTimepoints = size (all2all, 5);
bins  =  floor ( (nTimepoints/mf)- win_width/mf+1 );
trialN = size(all2all, 1);    
chanN = size(all2all, 3);


if meanInTime
    xM = zeros (trialN, bins,  chanN * length(f));
    yM = zeros (trialN, bins,  chanN * length(f));
else 
    xM = zeros (trialN, bins,  chanN * length(f) * win_width);
    yM = zeros (trialN, bins,  chanN * length(f) * win_width);
end
            
for timei = 1:bins 
    timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
    if meanInTime
        x = mean(all2all(:, 1,:,f, timeBins), 5, 'omitnan');
        x = reshape (x, [trialN, chanN * length(f)]);
        y = mean(all2all(:, 2,:,f,timeBins), 5, 'omitnan');
        y = reshape (y, [trialN, chanN * length(f)]);
    else
        x = all2all(:, 1,:,f,timeBins);
        x = reshape (x, [trialN, chanN * length(f)* win_width]);

        y = all2all(:, 2,:,f,timeBins);
        y = reshape (y, [trialN, chanN * length(f)* win_width]);
    end

    xM(:, timei, :) =  x;
    yM(:, timei, :) =  y;
end

fprintf('\n');

rsaZ = zeros (trialN, bins, bins);
parfor triali = 1:trialN
    mX= squeeze(xM(triali,:,:));
    mY= squeeze(yM(triali,:,:));
    r = corr (mX', mY','Type', 's'); 
    rsaZ(triali, :, :) = atanh(r);
end

rsaZ(isinf(rsaZ)) = nan;

end



 
 
 
 

