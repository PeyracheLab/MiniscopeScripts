% Pipeline for processing the miniscope data
% USAGE: MasterPreProcessing_Intan(folder_name)
% where folder_name is the base name for everything 

function Batch_CNMFE


fbasename = 'A6509';

basepath = '/media/guillaume/Elements/A6500/A6509';

file = '/home/guillaume/PSBImaging/python/datasets_A6509.csv';

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
        mergename = file_names{ii}
        
        %% Normcorre
        % takes the raw avi file
        avifile = fullfile(pwd, [mergename '_raw.avi']);
        Process_normcorre(avifile, mergename);
 
        %% CNMF-E
        neuron = Process_cnmfe(fullfile(pwd, [mergename '.h5']));

        % Exporting
        %save as csv        neuron = Process_cnmfe(fullfile(pwd, [mergename '.h5']));

        csvwrite(fullfile(pwd, [mergename '_C.csv']), neuron.C);
        csvwrite(fullfile(pwd, [mergename '_A.csv']), neuron.A);
        csvwrite(fullfile(pwd, [mergename '_C_raw.csv']), neuron.C_raw);
        csvwrite(fullfile(pwd, [mergename '_S.csv']), neuron.S);
        neuron.save_neurons();
        neuron.save_results([mergename '_cnmfe.mat']);
            
        %pwd
        cd('~')
        %[basepath, ~] = fileparts(dirName);
        cd(basepath);
        
        
        
    end
end


    

end
