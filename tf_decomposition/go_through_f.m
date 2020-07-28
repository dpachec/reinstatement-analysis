
function [] = go_through_f(win_size, step_size, onebyone)
folders = dir(); dirs = find(vertcat(folders.isdir));
folders = folders(dirs);


for foldi = 3:length(folders) %start at 3 cause 1 and 2 are . and ...
    
    direct = folders(foldi);
    

    cd (direct.name)
    sublist = dir('*all2all.mat');
    sublist = {sublist.name};


    [files] = load_files ('*all2all.mat');
    calculate_rsa(files, win_size, step_size, onebyone);

    % put all2all files in folders
    mkdir(direct.name);
    for filei=1:length(sublist)
        str = sublist{filei};
        movefile(str, direct.name)   
    end


    clearvars -except folders win_size step_size onebyone

    cd .. % goes up one directory

end   

toc
