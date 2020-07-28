
function [] = create_folders(fileNames)
sublist = dir(fileNames); sublist = {sublist.name};
fname_tmp = 'F0';
for filei=1:length(sublist)
   str = sublist{filei};
    fname_before = strsplit(str, '_');
    fname = fname_before{2};
    if strcmp(fname, 'all2alls.mat') fname = 'all2alls'; end
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

