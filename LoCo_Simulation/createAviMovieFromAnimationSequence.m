% createAviMovieFromAnimationSequence.m
%
% This script was copied and modified from the online Mathworks example for
% creating an .avi animation movie from a numbered sequence of .jpg images:
%
%    http://www.mathworks.com/help/matlab/examples/convert-between-image-sequences-and-video.html
%
% Marc Compere, comperem@gmail.com
% created : 11 Jan 2016
% modified: 17 Jan 2016

workingDir = 'anim_sequences'; % set this to wherever you want
                               % You might want to turn off Dropbox or Box
                               % or iCloud or Google Drive
                               % file synchronization services while making
                               % the animation. The simulation will go faster
                               % if it doesn't also have to sync all the
                               % new image files.

% -------------------------------------
% ------  Find Image File Names  ------
% -------------------------------------
% Find all the JPEG file names in the images folder. Convert the set of image names to a cell array.
disp(' ')
str=sprintf('using *all* images discovered in folder [%s]',workingDir); disp(str)
disp(' ')
imageNames = dir(fullfile(workingDir,'*.jpg'));
imageNames = {imageNames.name}';
str=sprintf('discovered [%i] image files for the animation',length(imageNames)); disp(str)

% figure out the name of the movie from the first filename
% If we've got a sequence of images named like this:
%   imageNames = 
%       '2016-01-11_17_06_55_img_000001.jpg'
%       '2016-01-11_17_06_55_img_000002.jpg'
%       '2016-01-11_17_06_55_img_000003.jpg'
%       '2016-01-11_17_06_55_img_000004.jpg'
%        ... (and so on)
% Then extract all occurrences of the '_' character and make the AVI
% filename everything up to the last underbar character.
filenameStr=imageNames{1};
k = strfind(filenameStr, '_');
AviFileName_prefix=filenameStr(1:(max(k)-1));
AviFileName = strcat(AviFileName_prefix,'.avi');
str=sprintf('using animation filename: [%s]',AviFileName); disp(str)

disp(' ')
pause(4)
disp(' ')

% --------------------------------------------------------
% ------  Create New Video with the Image Sequence  ------
% --------------------------------------------------------
% Construct a VideoWriter object, which creates a Motion-JPEG AVI file by default.
disp('creating a new VideoWriter() object...');
outputVideo = VideoWriter(fullfile(workingDir,AviFileName));
outputVideo.FrameRate = anim_fps;
open(outputVideo)
disp('done.');


% Loop through the image sequence, load each image, and then write it to the video.
disp('looping through each image in the sequence')
for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,imageNames{ii}));
   writeVideo(outputVideo,img)
   str=sprintf('done processing image [%i] of [%i], [%s]',ii,length(imageNames),imageNames{ii}); disp(str)
end
disp('done.');


% Finalize the video file.
disp('closing the outputVideo object')
close(outputVideo)
disp('done.');



disp('AVI animation file creation complete')
str=sprintf('you should now have a new .avi movie file in [%s] named [%s]',workingDir,AviFileName); disp(str)
% you can delete the .jpg sequence files once you've got the .avi.


% don't forget to turn your file sharing services back on like Dropbox or
% Google Drive
















