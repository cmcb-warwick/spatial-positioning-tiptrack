function [FileList] = loadData(loadDir, varargin)

% loadDir     = cellstr(['/Users/factory1/Desktop/Ben/Data'])
%  desc        = '*.tif';
%  dept        = 8;

n = 1;
k = 1;

% sets inputs, if no inputs provided, dept of intergation of files is
% parent and subfolders. and looks for tif files.

loadDir = cellstr(loadDir);

if length(varargin) > 0
    dept = varargin{1};
    if length(varargin) > 1
        desc = varargin{2};
    else
        desc = '/*.tif';
    end
else
    dept = 2;
end

% loads initial data;
LoadDir = loadDir;
for n = k:dept-1;
    for i = 1:length(loadDir)
        DirData              = dir(loadDir{i});
        DirIndex             = [DirData.isdir];
        subDir               = {DirData(DirIndex).name}';
        ValIndex(DirIndex)   = ~ismember(subDir,{'.','..','.DS_Store'});
        DirList              = {DirData(ValIndex).name}';
        for j = 1:length(DirList)
            DirList(j) = strcat(loadDir{i},'/',DirList(j));
        end
    LoadDir              = [LoadDir;DirList];
    clear ValIndex;
    end
    loadDir = LoadDir(k+1:length(LoadDir));
    k = length(LoadDir);
end
FileList = [];

for i = 1:length(LoadDir)
    files    = dir([LoadDir{i} desc]);
    files    = {files.name}';
    files    = sortrows(files);
    for j = 1:length(files)
        FileList = [FileList; strcat(LoadDir{i}, '/', files(j))];
    end
    
end



