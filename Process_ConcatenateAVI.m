% Pipeline for processing the miniscope data
% USAGE: MasterPreProcessing_Intan(folder_name)
% where folder_name is the base name for everything 

function Process_ConcatenateAVI(folder_name, spatial_downsampling)

%% Concatenate AVI files
%Process_ConcatenateAVI;

aviFiles = dir([folder_name filesep '*.avi']);
numFiles = length(aviFiles);

vidNum = [];
frameNum = [];
numFrames = 0;

%generate a vidObj for each video file. Also calculate total frames
for i=1:numFiles
    name = [folder_name filesep num2str(i-1) '.avi'];
    vidObj{i} = VideoReader(name);
    vidNum = [vidNum i*ones(1,vidObj{i}.NumberOfFrames)];
    frameNum = [frameNum 1:vidObj{i}.NumberOfFrames];
    numFrames = numFrames + vidObj{i}.NumberOfFrames;        
end
height = vidObj{1}.Height;
width = vidObj{1}.Width;

% go one folder up
idcs = strfind(folder_name, filesep);
newDir = folder_name(1:idcs(end)-1);
writerObj = VideoWriter([newDir filesep 'msvideo.avi'],'Grayscale AVI');
open(writerObj);

for video_i = 1:numFiles;
    name = [vidObj{1, video_i}.Path filesep vidObj{1, video_i}.Name];
    Yf = read_file(name);
    Yf = downsample_data(Yf,'space',1,spatial_downsampling,1);
    writeVideo(writerObj,uint8(Yf));
end

close(writerObj);



end