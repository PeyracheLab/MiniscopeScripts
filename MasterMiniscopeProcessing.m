% Pipeline for processing the miniscope data
% USAGE: MasterPreProcessing_Intan(folder_name)
% where folder_name is the base name for everything 

function MasterMiniscopeProcessing(fbasename,varargin)

%clear;
%fbasename = 'A0634';
%fbasename = 'Test'
%% Concatenate AVI files
spatial_downsampling = 2;

dirName = pwd;
avifolder = [dirName filesep 'Miniscope'];
Process_ConcatenateAVI(avifolder, spatial_downsampling);

%% CLEAN FOLDERS
% go one folder up and rename
idcs = strfind(dirName, filesep);
newDir = dirName(1:idcs(end)-1);
cd(newDir);
oldfoldername = dirName(idcs(end)+1:end);
tmp = erase(oldfoldername, '_');
newfoldername = [fbasename '-' tmp(3:end)];
mergename = [fbasename '-' tmp(3:end)];
movefile([newDir filesep oldfoldername], [newDir filesep newfoldername]);
cd([newDir filesep newfoldername]);

% rename analogin and tracking
files = dir('Take*.csv');
if exist(files.name, 'file')
    fname = files.name;
    targetFile = fullfile(pwd, [mergename '_0.csv']);
    movefile(fname, targetFile);
end

files = dir('Take*.tak');
if exist(files.name, 'file')
    fname = files.name;
    targetFile = fullfile(pwd, [mergename '_0.tak']);
    movefile(fname, targetFile);
end


% Copying analogin file if any exists 
intanfolder = dir([fbasename '_*']);
analoginfile = [intanfolder.name filesep 'analogin.dat'];
targetFile = fullfile(pwd, [mergename '_0_analogin.dat']);
movefile(analoginfile, targetFile);

% rename raw avi file
avifile = dir('*.avi');
movefile(avifile.name, fullfile(pwd, [mergename '_raw.avi']));

% rename json file
jsonfile = dir('*.json');
movefile(jsonfile.name, fullfile(pwd, [mergename '.json']));

% move miniscope timestamps in case
tstamps = dir('Miniscope/*.csv');
movefile(fullfile(tstamps.folder, tstamps.name), fullfile(pwd, [mergename '_ms_ts.csv']));

%cleaning
system(['rm -r ' fullfile(pwd, intanfolder.name)]);
system(['rm -r ' fullfile(pwd, 'experiment')]);

%% Normcorre
% takes the raw avi file
avifile = fullfile(pwd, [mergename '_raw.avi']);
Process_normcorre(avifile, mergename);

% Yf = read_file(fullfile(pwd, [mergename '.h5']));

%% CNMF-E
neuron = Process_cnmfe(fullfile(pwd, [mergename '.h5']));


%% Exporting
% save as csv
%csvwrite('A0634-201120_A.csv', neuron.A);
csvwrite(fullfile(pwd, [mergename '_C.csv']), neuron.C);
csvwrite(fullfile(pwd, [mergename '_A.csv']), neuron.A);
neuron.save_neurons();

end
