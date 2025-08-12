BASE_DIRECTORY = "..";

addpath(BASE_DIRECTORY);

DATE = '2025-08-08';
RUNINDEX = '21';

resFolder = fullfile(BASE_DIRECTORY, DATE, RUNINDEX, "RES");

% Get list of all subfolders in 'res'
subfolders = dir(fullfile(resFolder, "/RES_*"));

% Loop over each subfolder
for i = 1:length(subfolders)
    runName = subfolders(i).name;

    subfolderPath = fullfile(resFolder, runName);

    chunks = dir(fullfile(subfolderPath, "RES_*.mat"));
    chunkCount = length(chunks);
    
    invalidChunkCount = 0;
    for j = 1:chunkCount
        chunkName = chunks(j).name;

        matFilePath = fullfile(subfolderPath, chunkName);

        loadData = load(matFilePath);

        eta_max = max(max(max(abs(loadData.p.eta_data))));
        if eta_max > 0.1 || isnan(eta_max) || isinf(eta_max)
            invalidChunkCount = invalidChunkCount + 1;
            % break
        end
    end

    rating = "OK";
    if (invalidChunkCount == chunkCount)
        rating = "FAIL";
    elseif invalidChunkCount > 0
        rating = "MIXED";
    end

    fprintf("%s has rating %s.\n", runName, rating);
end