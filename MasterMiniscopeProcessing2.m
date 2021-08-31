cd % Pipeline for processing the miniscope data
% USAGE: MasterPreProcessing_Intan(folder_name)
% where folder_name is the base name for everything 

function MasterMiniscopeProcessing2(fbasename,varargin)


%% Normcorre
% takes the raw avi file
[~,mergename] = fileparts(pwd);

avifile = fullfile(pwd, [mergename '_raw.avi']);
Process_normcorre(avifile, mergename);

Y = read_file(fullfile(pwd, [mergename '.h5']));

%% CNMF-E
neuron = Process_cnmfe(fullfile(pwd, [mergename '.h5']));


%% Exporting
% save as csv
csvwrite(fullfile(pwd, [mergename '_C.csv']), neuron.C);
csvwrite(fullfile(pwd, [mergename '_A.csv']), neuron.A);
csvwrite(fullfile(pwd, [mergename '_C_raw.csv']), neuron.C_raw);
csvwrite(fullfile(pwd, [mergename '_S.csv']), neuron.S);
neuron.save_neurons();
neuron.save_results([mergename '_cnmfe.mat']);
end
