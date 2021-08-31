% Pipeline for processing the miniscope data
% USAGE: MasterPreProcessing_Intan(folder_name)
% where folder_name is the base name for everything 

function Batch_CNMFE


fbasename = 'A0634';

basepath = '/media/guillaume/Elements/A0600/A0634';

file = '/home/guillaume/PSBImaging/python/datasets_A0634.csv';

%%
fid = fopen(file);
for l=1:6
    tline = fgetl(fid);    
end
tline = fgetl(fid);
file_names = [];
while ischar(tline)
    if ~strcmp(tline(1), '#') & ~strcmp(tline(1), ',')    
        tmp = split(tline, ',');
        ymd = split(tmp{1}, '/');
        ymd = join(ymd, '');
        ymd = ymd{1};        
        file_names{end+1} = [fbasename '-' ymd(3:end)];
    end
    tline = fgetl(fid);
end


[n_folders,~] = size(file_names');

for ii=1:n_folders
    path = [basepath '/' file_names{ii}];
    if isfolder(path)
        cd(path)
        mergename = file_names{ii};

        %% Normcorre
        % takes the raw avi file
        avifile = fullfile(pwd, [mergename '_raw.avi']);
        Process_normcorre(avifile, mergename);
    
    end
end


    

end
