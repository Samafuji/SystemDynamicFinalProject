clc; clear; close all;

% User-defined parameters
n = 31; % Number of segments
num_parts = 5; % For symmetric partition
E = 27.8e9; % Young's modulus (N/m^2)
D = 2400; % Density (kg/m^3)
W = 35; % Width (m)
T = 9; % Thickness (m)
L = 1700; % Total length of the float bridge (m)

alpha = 0.003;
Kground = 1e10;
Cground = alpha * Kground;
Kanchor = 0.7e8;
M_box = 11000000; % Mass of box for pier
waveFrequency = 0.1; % Wave frequency (Hz)
windFrequency = 4.5963; % Wind frequency (Hz)
simTime = 1000; % Total simulation time (s)
dt = 0.1;
tspan = 0:dt:simTime;

% Build bridge configuration (0 or 1 for pier presence)
result = symmetric_partition(n, num_parts)

wight = D*W*T*L+M_box*num_parts

% Horizontal system matrices
[A_horizontal, B_horizontal, C_out_horizontal, D_horizontal, ~, ~, ~, m] = ...
    calculate_system_matrices(n, result, E, D, W, T, L, alpha, Kground, Cground, Kanchor, M_box, 1);

% Vertical system matrices
[A_vertical, B_vertical, C_out_vertical, D_vertical, ~, ~, ~, ~] = ...
    calculate_system_matrices(n, result, E, D, W, T, L, alpha, Kground, Cground, Kanchor, M_box, 0);

% Calculate external forces
[Fx, Fy] = calculate_forces(n, result, windFrequency, waveFrequency, L/n, T, W, D, tspan, m, M_box);

Fx(:,:) = 0;

% Initial conditions
x0 = zeros(2 * n, 1);

% Run horizontal simulation
sys_horizontal = ss(A_horizontal, B_horizontal, C_out_horizontal, D_horizontal);
[y_horizontal, ~, ~] = lsim(sys_horizontal, Fx, tspan, x0);

% Run vertical simulation
sys_vertical = ss(A_vertical, B_vertical, C_out_vertical, D_vertical);
[y_vertical, ~, ~] = lsim(sys_vertical, Fy, tspan, x0);

% Unified 3D Visualization
visualize_results_unified_3d(tspan, y_horizontal, y_vertical, n, L, result);




% -------------------------------------------------------------------------
% Calculate System Matrices
% -------------------------------------------------------------------------
function [A, B, C_out, D_M, M_Matrix, n, l, m] = calculate_system_matrices(n, result, E, D, W, T, L, alpha, Kground, Cground, Kanchor, M_box, yoko)
    l = L / n;      % segment length
    m = l * T * W * D; % mass per block
    
    % Second moment of inertia (assume rectangular cross-section)
    if yoko == 1
        % Horizontal axis (yoko == 1)
        Ia = (W^3 * T) / 12; 
    else
        % Vertical axis (yoko == 0)
        Ia = (W * T^3) / 12; 
    end
    
    % Baseline stiffness
    Kr = 8 * E * Ia * l / l^4;  
    Cr = alpha * Kr;

    M_list = ones(1, n) * m; 
    C_list = ones(1, n) * Cr; 
    K_list = ones(1, n) * Kr; 

    CWater = 0; 
    KWater = 0;
    
    % Mass Matrix
    M_Matrix = diag(M_list);
    for i = 1:n
        if result(i) == 1
            M_Matrix(i, i) = M_Matrix(i, i) + M_box;
        end
    end

    % Damping Matrix
    C_Matrix = zeros(n);
    for i = 1:n
        if i > 1
            C_Matrix(i, i-1) = -C_list(i);
            C_Matrix(i-1, i) = -C_list(i);
        end
        if i == 1 || i == n
            C_Matrix(i, i) = C_list(i) + Cground;
        elseif result(i) == 1
            C_Matrix(i, i) = C_list(i) + C_list(i+1) + CWater;
        else
            C_Matrix(i, i) = C_list(i) + C_list(i+1);
        end
    end

    % Stiffness Matrix
    K_Matrix = zeros(n);
    for i = 1:n
        if i > 1
            K_Matrix(i, i-1) = -K_list(i);
            K_Matrix(i-1, i) = -K_list(i);
        end
        if i == 1 || i == n
            K_Matrix(i, i) = K_list(i) + Kground;
        elseif result(i) == 1
            K_Matrix(i, i) = K_list(i) + K_list(i+1) + KWater;
        else
            K_Matrix(i, i) = K_list(i) + K_list(i+1) + Kanchor;
        end
    end

    % State-space matrices
    A = [zeros(n), eye(n); -M_Matrix \ K_Matrix, -M_Matrix \ C_Matrix];
    B = [zeros(n); inv(M_Matrix)];
    C_out = eye(2 * n); 
    D_M = zeros(2 * n, n);
end


% -------------------------------------------------------------------------
% Calculate Input Forces (wind and wave)
% -------------------------------------------------------------------------
function [Fx, Fy] = calculate_forces(n, result, windFrequency, waveFrequency, l, T, W, D, tspan, m, M_box)
    % Wind force calculation
    % Generate Gaussian data
    g = 9.81;
    mu = 16; sigma = 1;
    data = (mu + sigma * randn(length(tspan), 1));
    u = zeros(length(tspan), n);
    for i = 1:n
        u(:, i) = data*sin(i*pi/n);
    end

    Cd = 0.7;
    Area = l * T; 
    low = 1.293;
    Fwind = 1/2 * low * Area * u.^2 * Cd;

    % Wave forces
    [FwaveX_basic, FwaveY_basic] = calculate_wave_forces(waveFrequency, 5, 10, tspan);

    ux = zeros(length(tspan), n);
    uy = zeros(length(tspan), n);
    for i = 1:n
        if result(i) == 1
            ux(:, i) = FwaveX_basic;
            uy(:, i) = FwaveY_basic - m * g - M_box * g;
        else
            uy(:, i) = -m*g;
        end
    end

    Fx = Fwind + ux;
    Fy = uy;
end

% -------------------------------------------------------------------------
% System Properties
% -------------------------------------------------------------------------
function [eigenvalues, natural_frequencies, damping_ratios] = system_properties(A)
    eigenvalues = eig(A);
    natural_frequencies = abs(eigenvalues) / (2 * pi);
    damping_ratios = -real(eigenvalues) ./ abs(eigenvalues);
end

% -------------------------------------------------------------------------
% Visualization
% -------------------------------------------------------------------------
function visualize_results(t, y, n, L, result)
    bridge_x = linspace(0, L, n);
    bridge_y = zeros(1, n);
    zMax = max(max(y(:, 1:n))) * 1.01;

    figure;
    plot(t, y(:, 1:n));
    xlabel('Time (s)');
    ylabel('Displacement');
    title('Displacement Response at Each Node');
    grid on;
    legend;

    figure;
    plot(t, y(:, n+1:end));
    xlabel('Time (s)');
    ylabel('Velocity');
    title('Velocity Response at Each Node');
    grid on;
    legend;

    if 1
        figure(99);
        for i = 1:length(t)
            bridge_z = y(i, 1:n);
            plot3(bridge_x, bridge_y, bridge_z, '-o', 'LineWidth', 2);
            hold on;
            scatter3(bridge_x(result == 1), bridge_y(result == 1), bridge_z(result == 1), 50, 'r', 'filled');
            hold off;
            axis([0 L -2 2 -zMax zMax]);
            xlabel('Bridge Length (m)');
            ylabel('Width (m)');
            zlabel('Displacement (m)');
            title(sprintf('Time: %.2f seconds', t(i)));
            grid on;
            pause(0.001);
        end
    end
end

function visualize_results_unified_3d(t, y_horizontal, y_vertical, n, L, result)
    % Calculate x positions along the bridge
    bridge_x = linspace(0, L, n);

    % Calculate maximum displacements for both axes
    yMax = max(abs(y_horizontal(:, 1:n)), [], 'all') * 1.1; % Max horizontal displacement
    zMax = max(abs(y_vertical(:, 1:n)), [], 'all') * 1.1;   % Max vertical displacement

    % Ensure positive values for axis limits
    if yMax == 0, yMax = .1; end
    if zMax == 0, zMax = .1; end

    % Create figure for unified 3D animation
    figure(99);
    for i = 1:length(t)
        % Extract displacements at the current time step
        bridge_y_horizontal = y_horizontal(i, 1:n); % Horizontal displacements
        bridge_z_vertical = y_vertical(i, 1:n);     % Vertical displacements

        % Plot unified 3D bridge animation
        plot3(bridge_x, bridge_y_horizontal, bridge_z_vertical, '-o', 'LineWidth', 2);
        hold on;

        % Highlight piers
        scatter3(bridge_x(result == 1), bridge_y_horizontal(result == 1), ...
                 bridge_z_vertical(result == 1), 50, 'r', 'filled'); % Piers
        hold off;

        % Axis settings
        axis([0 L -yMax yMax -zMax zMax]); % Ensure valid ranges for axis limits
        xlabel('Bridge Length (m)');
        ylabel('Horizontal Displacement (m)');
        zlabel('Vertical Displacement (m)');
        title(sprintf('Unified 3D Displacement Animation\nTime: %.2f seconds', t(i)));
        grid on;

        pause(0.001); % Pause for animation effect
    end
end
