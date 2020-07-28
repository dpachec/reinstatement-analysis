
function [files] = load_files(files2load)

sublist = dir(files2load);
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subno=1:length(sublist)
    disp (['File nº ' num2str(subno)]);
    load(sublist{subno});
    s_all{subno} = all2all; 
    %size(all2all)
    if exist ('chanNames')
        chanNames_all{subno} = chanNames;
    else
        chanNames_all = [];
    end
end

if exist('s_all')
    files               = [];
    files.s_all         = s_all; 
    files.chanNames_all = chanNames_all; 
    files.sublist       = sublist;
else
   files = []; 
end


