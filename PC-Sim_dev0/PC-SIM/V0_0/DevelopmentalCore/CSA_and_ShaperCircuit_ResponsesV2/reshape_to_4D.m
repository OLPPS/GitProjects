function [V4, Nx, Ny, Nz] = reshape_to_4D(data)
% data is MxN: [x, y, z, v1, v2, ..., vT]

    % Extract coordinates
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);

    [x1,x2] = maxmin(x)
    [y1,y2] = maxmin(y)
    [z1,z2] = maxmin(z)

    % Extract values
    V = data(:,4:end);   % M x T
    T = size(V,2);

    % Unique coordinate grids
    ux = unique(x);
    uy = unique(y);
    uz = unique(z);

    Nx = numel(ux);
    Ny = numel(uy);
    Nz = numel(uz);

    % Preallocate output
    V4 = nan(T, Nx, Ny, Nz);

    % Build index maps
    [~, ix] = ismember(x, ux);
    [~, iy] = ismember(y, uy);
    [~, iz] = ismember(z, uz);

    % Fill the 4D array
    for row = 1:length(x)
        V4(:, ix(row), iy(row), iz(row)) = V(row, :).';
    end
end
