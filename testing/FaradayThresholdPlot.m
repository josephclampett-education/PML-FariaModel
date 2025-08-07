BASE_DIRECTORY = "..";

addpath(BASE_DIRECTORY);

% Get list of all .mat files in the folder
files = dir(fullfile("../threshold_cache", '*.mat'));

% Initialize arrays to store h0, h1, and R
h0_values = [];
h1_values = [];
R_values = [];
GamF_values = [];

for i = 1:length(files)
    % Get filename without extension
    [folder, name, ext] = fileparts(files(i).name);
    
    % Split the name by underscore
    parts = split(name, '_');
    
    if length(parts) == 3
        % Convert each part to number
        h0 = str2double(parts{1});
        h1 = str2double(parts{2});
        R  = str2double(parts{3});
        
        % Store the values
        if (h0 ~= h1)
            h0_values(end+1) = h0 * 10^3;
            h1_values(end+1) = h1 * 10^3;
            R_values(end+1)  = R;
    
            GamF_values(end+1) = load(fullfile("../threshold_cache", files(i).name)).GamF;
        end
    else
        warning('Filename "%s" does not match expected format.', name);
    end
end


% Create grid (only if your h1 and R form a complete grid)
[h1_unique, ~, h1_idx] = unique(h1_values);
[R_unique, ~, R_idx] = unique(R_values);

% Initialize grid for GamF values
GamF_grid = NaN(length(R_unique), length(h1_unique));

% Fill in the GamF grid
for i = 1:length(GamF_values)
    GamF_grid(R_idx(i), h1_idx(i)) = GamF_values(i);
end

% Create meshgrid
[H1, Rgrid] = meshgrid(h1_unique, R_unique);

% Plot surface
figure;
surf(H1, Rgrid, GamF_grid);
hold on
plot3(h1_values, R_values, GamF_values,'o')
xlabel('h1');
ylabel('R');
zlabel('GamF');
title('Surface plot of GamF vs h1 and R');
colorbar;
shading interp;  % Optional: smoothens the surface

saveas(gcf, "ThresholdMap");
