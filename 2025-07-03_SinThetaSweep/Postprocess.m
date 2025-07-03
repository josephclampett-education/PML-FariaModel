fileSearch = dir('*.mat');
fileCount = length(fileSearch);

sinThetas = zeros(fileCount, 1);
meanVelocities = zeros(fileCount, 1);
for i = 1:fileCount
    filePath = fullfile(fileSearch(i).folder, fileSearch(i).name);
    load(filePath);

    dx = diff(p.x_data);
    dy = diff(p.y_data);
    dt = diff(p.t_data);
    
    v_data = sqrt(dx.^2 + dy.^2) ./ dt * 4.75 * 40;
    meanVelocity = mean(v_data);
    
    sinThetas(i) = p.sin_theta;
    meanVelocities(i) = meanVelocity;
end
scatter(sinThetas, meanVelocities);
title('Droplet Speed (mm/s) vs. $\sin{\theta}$', 'Interpreter', 'latex')
xlabel('$\sin{\theta}$','Interpreter','latex')
ylabel('$v$ (mm/s)','Interpreter','latex')