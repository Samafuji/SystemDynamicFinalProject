clc; clear; close all;

% User-defined parameters
n = 501; % Number of segments
num_parts = 5; % For symmetric partition
E = 27.8e9; % Young's modulus (N/m^2)
D = 2400; % Density (kg/m^3)
W = 35; % Width (m)
T = 9; % Thickness (m)
L = 1700; % Total length of the float bridge (m)

alpha = 0.3;
Kground = 1e11; 
Cground = alpha * Kground;
Kanchor = 0.7e8;
M_box = 18984932.7; % Mass of box for pier
waveFrequency = 0.3; % Wave frequency (Hz)
windFrequency = 4.5963; % Wind frequency (Hz)
simTime = 1000; % Total simulation time (s)
dt = 0.1;
tspan = 0:dt:simTime;

% Build bridge configuration (0 or 1 for pier presence)
result = symmetric_partition(n, num_parts)

wight = D*W*T*L+M_box*num_parts
% wight = D*W*T*L;



% Horizontal system matrices
[A_horizontal, B_horizontal, C_out_horizontal, D_horizontal, ~, ~, ~, m] = ...
    calculate_system_matrices(n, result, E, D, W, T, L, alpha, Kground, Cground, Kanchor, M_box, 1);

% Vertical system matrices
[A_vertical, B_vertical, C_out_vertical, D_vertical, ~, ~, ~, ~] = ...
    calculate_system_matrices(n, result, E, D, W, T, L, alpha, Kground, Cground, Kanchor, M_box, 0);

% Calculate external forces
[Fx, Fy] = calculate_forces(n, result, windFrequency, waveFrequency, L/n, T, W, D, tspan, m, M_box,wight);

% Fy(:,:) = 0;

% Initial conditions
x0 = zeros(2 * n, 1);

% Run horizontal simulation
sys_horizontal = ss(A_horizontal, B_horizontal, C_out_horizontal, D_horizontal);
[y_horizontal, ~, ~] = lsim(sys_horizontal, Fx, tspan, x0);


% Run vertical simulation
sys_vertical = ss(A_vertical, B_vertical, C_out_vertical, D_vertical);
[y_vertical, ~, ~] = lsim(sys_vertical, Fy, tspan, x0);

eigenvalues = eig(A_vertical);
natural_frequencies = abs(eigenvalues) / (2 * pi);
damping_ratios = -real(eigenvalues) ./ abs(eigenvalues);

% Unified 2D Visualization
visualize_results(tspan, y_horizontal, y_vertical, n, L, result);

% Unified 3D Visualization
visualize_results_unified_3d(tspan, y_horizontal, y_vertical, n, L, result);

% bode plot
nodes_to_plot = [1, 101, 201, 251];
bodePlot(sys_horizontal, sys_vertical, nodes_to_plot);

impulseResponsePlot(sys_horizontal, sys_vertical, nodes_to_plot);

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
    Kr = E * Ia * l / l^4;  
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
        elseif yoko == 0
            K_Matrix(i, i) = K_list(i) + K_list(i+1) + Kanchor;
        else
            K_Matrix(i, i) = K_list(i) + K_list(i+1);
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
function [Fx, Fy] = calculate_forces(n, result, windFrequency, waveFrequency, l, T, W, D, tspan, m, M_box,wight)
    % Wind force calculation
    % Generate Gaussian data
    g = 9.81;
    mu = 25; sigma = 1;
    data = (mu + sigma * randn(length(tspan), 1));
    u = zeros(length(tspan), n);
    for i = 1:n
        u(:, i) = data;
    end

    Cd = 0.7;
    Area = l * T; 
    low = 1.293;
    Fwind = 1/2 * low * Area * u.^2 * Cd;

    % Wave forces
    [FwaveX_basic, FwaveY_basic] = calculate_wave_forces(waveFrequency, 5, 10, tspan,wight);

    ux = zeros(length(tspan), n);
    uy = zeros(length(tspan), n);
    Fground = 6334189235;
    for i = 1:n
        if result(i) == 1
            ux(:, i) = FwaveX_basic;
            uy(:, i) = FwaveY_basic - m * g - M_box * g;
        elseif i == 1 || i == n
            uy(:, i) = -m*g + Fground;
        else
            uy(:, i) = -m*g;
        end
    end

    Fx = Fwind + ux;
    Fy = uy;
end

% -------------------------------------------------------------------------
% Visualization bodePlot
% -------------------------------------------------------------------------
function bodePlot(sys_horizontal, sys_vertical, list)
    for i = 1:length(list)
        node_idx = list(i);
        
        % Horizontal system Bode plot
        sub_sys_horizontal = sys_horizontal(node_idx, node_idx);
        figure;
        bode(sub_sys_horizontal);
        grid on;
        title(sprintf('Bode Plot (Horizontal) for Node %d', node_idx));

        % Vertical system Bode plot
        sub_sys_vertical = sys_vertical(node_idx, node_idx);
        figure;
        bode(sub_sys_vertical);
        grid on;
        title(sprintf('Bode Plot (Vertical) for Node %d', node_idx));
    end
end

% -------------------------------------------------------------------------
% Visualization impulseResponsePlot
% -------------------------------------------------------------------------
function impulseResponsePlot(sys_horizontal, sys_vertical, list)
    for i = 1:length(list)
        node_idx = list(i);

        % Horizontal impulse response
        figure;
        impulse(sys_horizontal(node_idx, node_idx));
        title(sprintf('Impulse Response (Horizontal) for Nodes: %s', node_idx));
        xlabel('Time (s)');
        ylabel('Amplitude (m)');
        grid on;
    
        % Vertical impulse response
        figure;
        impulse(sys_vertical(node_idx, node_idx));
        title(sprintf('Impulse Response (Vertical) for Nodes: %s', node_idx));
        xlabel('Time (s)');
        ylabel('Amplitude (m)');
        grid on;
    end
end

% -------------------------------------------------------------------------
% Visualization phasePortrait
% -------------------------------------------------------------------------
function phasePortrait(sys_horizontal, sys_vertical, nodes, tspan, x0)
    % Simulate the system response for initial conditions
    [y_horizontal, ~, ~] = lsim(sys_horizontal, zeros(length(tspan), size(sys_horizontal.B, 2)), tspan, x0);
    [y_vertical, ~, ~] = lsim(sys_vertical, zeros(length(tspan), size(sys_vertical.B, 2)), tspan, x0);

    % Plot phase portraits for horizontal system
    figure;
    for i = 1:length(nodes)
        node = nodes(i);
        subplot(ceil(sqrt(length(nodes))), ceil(sqrt(length(nodes))), i);
        plot(y_horizontal(:, node), y_horizontal(:, node + size(sys_horizontal.A, 1)/2));
        xlabel('Displacement (m)');
        ylabel('Velocity (m/s)');
        title(sprintf('Horizontal Phase Portrait for Node %d', node));
        grid on;
    end
    suptitle('Horizontal System Phase Portraits');

    % Plot phase portraits for vertical system
    figure;
    for i = 1:length(nodes)
        node = nodes(i);
        subplot(ceil(sqrt(length(nodes))), ceil(sqrt(length(nodes))), i);
        plot(y_vertical(:, node), y_vertical(:, node + size(sys_vertical.A, 1)/2));
        xlabel('Displacement (m)');
        ylabel('Velocity (m/s)');
        title(sprintf('Vertical Phase Portrait for Node %d', node));
        grid on;
    end
    suptitle('Vertical System Phase Portraits');
end
% -------------------------------------------------------------------------
% Visualization 2D
% -------------------------------------------------------------------------
function visualize_results(t, y_horizontal, y_vertical, n, L, result)
    % Define the node indices to plot
    nodes_to_plot = [1, 251, n, find(result == 1)]; % i == 1, i == 251, i == n, and piers
    nodes_to_plot = unique(nodes_to_plot); % Remove duplicates
    
    % Horizontal Displacement Plot
    figure;
    plot(t, y_horizontal(:, nodes_to_plot));
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    title('Horizontal Displacement Response (Selected Nodes)');
    grid on;
    legend(arrayfun(@(x) sprintf('Node %d', x), nodes_to_plot, 'UniformOutput', false));
    
    % Horizontal Velocity Plot
    figure;
    plot(t, y_horizontal(:, nodes_to_plot + n)); % Velocity components
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Horizontal Velocity Response (Selected Nodes)');
    grid on;
    legend(arrayfun(@(x) sprintf('Node %d', x), nodes_to_plot, 'UniformOutput', false));
    
    % Vertical Displacement Plot
    figure;
    plot(t, y_vertical(:, nodes_to_plot));
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    title('Vertical Displacement Response (Selected Nodes)');
    grid on;
    legend(arrayfun(@(x) sprintf('Node %d', x), nodes_to_plot, 'UniformOutput', false));
    
    % Vertical Velocity Plot
    figure;
    plot(t, y_vertical(:, nodes_to_plot + n)); % Velocity components
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Vertical Velocity Response (Selected Nodes)');
    grid on;
    legend(arrayfun(@(x) sprintf('Node %d', x), nodes_to_plot, 'UniformOutput', false));

end


% -------------------------------------------------------------------------
% Visualization 3D
% -------------------------------------------------------------------------
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
    for i = 1:length(t)
        figure(99);
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
        axis([0 L -zMax zMax -zMax zMax]); % Ensure valid ranges for axis limits
        % axis([0 L -yMax yMax -zMax zMax]); % Ensure valid ranges for axis limits
        xlabel('Bridge Length (m)');
        ylabel('Horizontal Displacement (m)');
        zlabel('Vertical Displacement (m)');
        title(sprintf('Unified 3D Displacement Animation\nTime: %.2f seconds', t(i)));
        grid on;

        pause(0.1); % Pause for animation effect
    end
end
