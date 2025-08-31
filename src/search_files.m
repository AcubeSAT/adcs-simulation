%% ======================================================================= %%
%  This script finds all files that contain a specific word or phrase.
%  Useful for knowing in which files a specific variable or function is used.
%  It is an auxiliary tool, not essential for running simulations.
% ======================================================================== %%


% Define the root directory
rootDir = '..';

% Get the filename of the current script file
scriptFilename = mfilename;

% Define the search phrase
searchPhrase = 'Q_selection';

% Recursively search for files in all subdirectories
files = dir(fullfile(rootDir, '**', '*.m'));

% Initialize a cell array to store results
matchingFiles = {};

% Loop through each file and search for the phrase
for i = 1:numel(files)
    filePath = fullfile(files(i).folder, files(i).name);
    
    [~, fileName, ~] = fileparts(filePath);
    if ~strcmp(fileName, scriptFilename)
        % Read the contents of the file
        fid = fopen(filePath, 'r');
        fileContent = fread(fid, '*char')';
        fclose(fid);
        
        % Search for the phrase in the file content
        if contains(fileContent, searchPhrase)
            matchingFiles{end+1} = filePath;
        end
    
    end

end

% Display the matching files
disp(['Matching files for search phrase: "', searchPhrase, '":']);

for i = 1:numel(matchingFiles)
    disp(matchingFiles{i});
end