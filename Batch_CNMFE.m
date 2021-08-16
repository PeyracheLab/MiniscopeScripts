% Pipeline for processing the miniscope data
% USAGE: MasterPreProcessing_Intan(folder_name)
% where folder_name is the base name for everything 

function Batch_CNMFE


fbasename = 'A0634';
file = '/media/Hypernova/Guillaume/A0600/A0634/datasets_A0634.csv';

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

basepath = fileparts(file);

[n_folders,~] = size(file_names');

for ii=1:n_folders
    path = [basepath '/' file_names{ii}];
    if isfolder(path)
        cd(path)
        mergename = file_names{ii};
%         %% Cleaning folder of previous cnmfe        
%          previous_init = [path filesep [mergename '_source_extraction']];
% 
%         filePattern = fullfile(previous_init, '*.fig'); % Change to whatever pattern you need.
%         theFiles = dir(filePattern);
%         for k = 1 : length(theFiles)
%           baseFileName = theFiles(k).name;
%           fullFileName = fullfile(myFolder, baseFileName);
%           fprintf(1, 'Now deleting %s\n', fullFileName);
%           delete(fullFileName);
%         end
%          
        %% Downsampling
        spatial_downsampling = 2;
        avifolder = [path filesep 'Miniscope'];
        Process_ConcatenateAVI(avifolder, spatial_downsampling);
        movefile('msvideo.avi', [mergename '_raw.avi']);
        
        %% Normcorre
        % takes the raw avi file
        avifile = fullfile(pwd, [mergename '_raw.avi']);
        Process_normcorre(avifile, mergename);

        % Yf = read_file(fullfile(pwd, [mergename '.h5']));
    
        
        %% CNMF-E
        try 
            neuron = Process_cnmfe(fullfile(pwd, [mergename '.h5']));

            % Exporting
            %save as csv
            csvwrite(fullfile(pwd, [mergename '_C.csv']), neuron.C);
            csvwrite(fullfile(pwd, [mergename '_A.csv']), neuron.A);
            csvwrite(fullfile(pwd, [mergename '_C_raw.csv']), neuron.C_raw);
            csvwrite(fullfile(pwd, [mergename '_S.csv']), neuron.S);
            neuron.save_neurons();
            neuron.save_results([mergename '_cnmfe.mat']);
        catch
            fprintf('error in cnmfe');
        end
            
            
        %pwd
        cd('~')
        %[basepath, ~] = fileparts(dirName);
        cd(basepath);

    end
end


    

end
