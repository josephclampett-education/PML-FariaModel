BASE_DIRECTORY = "../..";

addpath(BASE_DIRECTORY);

% ================================================================
% PARAMS
% ================================================================

% Saving
VAR_outputFolder = "RES";

%% ================================================================

folders = dir(fullfile(VAR_outputFolder, "RES_*"));
threadCount = length(folders);

outputData = [];

parfor i = 1:threadCount
    folderPath = fullfile(folders(i).folder, folders(i).name);

    % pathId is important for saving!
    [pathstr, name, ext] = fileparts(folderPath);
    pathId = fullfile(folderPath, strcat(name, ext)); % Add duplicate of innermost folder string

    fileSearchPath = fullfile(folderPath, "*.mat");
    files = dir(fileSearchPath);

    for j = 1:numel(files)
        filePath = fullfile(files(j).folder, files(j).name);
        loadP = load(filePath);
        p = loadP.p;

        x_data = p.x_data;
        y_data = p.x_data;

        output_filename = files(j).name;
        % generate a video for each res_i
        v = VideoWriter([output_filename, '.avi']);
        v.FrameRate = 40;
        open(v);

        % Create a figure without displaying it (for cluster use)
        fig = figure('Visible', 'off');
        axis tight manual

        for t = 1:M
            clf; % Clear the figure
            scatter(x_data(t, :), y_data(t, :), 36, 'filled');
            title(sprintf('Frame %d / %d', t, M));
            xlim(x_range);
            ylim(y_range);
            xlabel('X');
            ylabel('Y');
            drawnow;

            % Capture the frame
            frame = getframe(fig);
            writeVideo(v, frame);
        end

        close(v);
        close(fig); % Clean up

        fprintf('Video saved to %s\n', output_filename);

    end
    
end