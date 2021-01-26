function [ flist ] = get_file_list( data_dir, ext, search_str )
%GET_FILE_LIST read list of files in directory
% INPUT:
%       data_dir: directory to list; current directory will not change
%           after function call (default is pwd)
%       ext: extension of files to list, as string, without leading dot (defauls is '*')
% OUTPUT:
%       flist: cell array of files in given directory
%
% Kristof Bognar May 2018

%% define directory and extension if no input

changedir=true;

if nargin==0
    data_dir=pwd();     % use current directory
    ext='*';            % read all files
    changedir=false;    % don't try changing directories
    search_str='*';
elseif nargin==1
    ext='*';
    search_str='*';
elseif nargin==2
    search_str='*';
end


%% read list of files

% change directory if required
if changedir
    cur_dir=pwd();
    cd(data_dir);
end

% list files
tmp=dir([search_str '.' ext]);
flist={tmp.name};

% check if any files exist
if isempty(flist), return, end

% remove . and .. if extension is not specified
if strcmp(flist(1),'.')
    flist(1)=[];
end

if strcmp(flist(1),'..')
    flist(1)=[];
end

% remove temporary backup files (*.ext~)
ind_tmp = strfind(flist, '~');
ind_remove = find(not(cellfun('isempty', ind_tmp)));
flist(ind_remove)=[];

% return to original directory
if changedir
    cd(cur_dir);
end

end

