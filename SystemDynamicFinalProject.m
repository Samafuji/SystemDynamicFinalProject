clc; clear; close all;

% ----------------------------------------------------------------
% Build Bridge
% ----------------------------------------------------------------
n = 31;
num_parts = 5; 
result = symmetric_partition(n, num_parts);

% or
% result = [0 0 1 0 0 1 0 0 1 0 0];

% ----------------------------------------------------------------
% PROPERTIES
% ----------------------------------------------------------------
n = length(result);

%RC Concrete
E = 27.8*10^9;  %Young's modulus (K/m^2)
D = 2400;       %Density (Kg/m^3)
W = 35;        %width  m
T = 9;       %thick  m
L = 1700;        %Total length of the float bridge  m

l = L/n;        %length of the bridge per unit
m = 1*0.3*0.04*2400;  % mass of the bridge per unit
m = l * T * W * D; % per block
Ia = W*T^3/12;  %second moment of inertia
Kr = E * Ia /l^4; %stiffness coefficient K = E * Ia / l; ==> F = kr * theta
Kr = 8 * E * Ia * l /l^2; % x = -F_bar * l^4 / (8 *E * Ia) && F / l = F_bar
alpha = 0.003;
Cr = alpha*Kr; %damping coefficient

%  23m
% ----
% |  | 110m  11000 tons  height 8.5m
% ----



% Matrices
M_list = ones(1, n)* m; % Mass values
C_list = ones(1, n) * Cr; % Damping values
K_list = ones(1, n) * Kr; % Stiffness values

CWater = 0; Kground = 1e7; Cground = alpha * Kground; KWater = 0;
M_box = 11000000; % 11000 tons
K_anchor = 0;

% M_Matrix
M_Matrix = diag(M_list);
for i = 1:n
    if result(i) == 1
        M_Matrix(i,i) = M_Matrix(i,i) + M_box;
    end
end

% C_Matrix
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

% K_Matrix
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
        K_Matrix(i, i) = K_list(i) + K_list(i+1) + K_anchor;
    end
end

%% State-space representation
A = [zeros(n), eye(n); -M_Matrix \ K_Matrix, -M_Matrix \ C_Matrix];
B = [zeros(n); inv(M_Matrix)];
C_out = eye(2 * n); % Identity for full-state output
D = zeros(2 * n, n); % No direct transmission

% %% transfer funciton
% [b, a] = ss2tf(A,B,C_out,D,1);
% TF = tf(b,a)

% ----------------------------------------------------------------
% Input PROPERTIES
% ----------------------------------------------------------------

f = 4.5963; % Frequency of the sine wave (Hz)
omega = 2 * pi * f; % Angular frequency
tspan = 0:0.1:1000; % Time span
% Generate sinusoidal input at a specific node
u_node = round(n / 2); % Apply the sine wave at the middle node
u = zeros(length(tspan), n); % Initialize input matrix
for t = 1:length(tspan)
    u(t, u_node+1) = 30 * (sin(omega * tspan(t)) + 1); % Sinusoidal input at the middle node
end


% 生成高斯分布數據
mu = 16;          % 均值
sigma = 1;       % 標準差
data = (mu + sigma * randn(length(tspan), 1));
for i = 1:n
   u(:, i) = data*sin(i*pi/n);
end

v = u;
Cd = 0.7;
Area = l * T; % Area = length of bridge * width / n
low = 1.293;
Fwind = 1/2 * low * Area * v.^2 * Cd;

% Initial conditions
x0 = zeros(2 * n, 1);

sys = ss(A, B, C_out, D); % State-space system
[y, t, x] = lsim(sys, Fwind, tspan, x0);
% [y, t] = impulse(sys);

% Plot displacement response
figure;
plot(t, y(:, 1:n));
xlabel('Time (s)');
ylabel('Displacement');
title('Displacement Response at Each Node (Sine Wave Input)');
grid on;
legend;
set(gcf, 'Position', [80, 250, 800, 400])

% Plot velocity response
figure;
plot(t, y(:, n+1:end));
xlabel('Time (s)');
ylabel('Velocity');
title('Velocity Response at Each Node (Sine Wave Input)');
grid on;
legend;
set(gcf, 'Position', [1050, 250, 800, 400])


% Eigenvalues of the system
eigenvalues = eig(A);

% Natural frequencies (in Hz)
natural_frequencies = abs(eigenvalues) / (2 * pi); % Magnitude of eigenvalues converted to Hz

% Damping ratios (optional)
damping_ratios = -real(eigenvalues) ./ abs(eigenvalues);

% Display results
disp('Eigenvalues:');
disp(eigenvalues);

disp('Natural Frequencies (Hz):');
disp(natural_frequencies);

disp('Damping Ratios:');
disp(damping_ratios);

result

% ----------------------------------------------------------------
% Visualization: 3D Animation
% ----------------------------------------------------------------
bridge_x = linspace(0, L, n);
bridge_y = zeros(1, n);

figure(99);
for i = 1:length(t)
    % Update bridge displacement
    bridge_z = y(i, 1:n);
    
    % Plot 3D bridge
    plot3(bridge_x, bridge_y, bridge_z, '-o', 'LineWidth', 2);
    hold on;
    scatter3(bridge_x(result == 1), bridge_y(result == 1), bridge_z(result == 1), 50, 'r', 'filled'); % with pier
    hold off;
    
    % Axis settings
    axis([0 L -2 2 -0.1 0.1]);
    xlabel('Bridge Length (m)');
    ylabel('Width (m)');
    zlabel('Displacement (m)');
    title(sprintf('Time: %.2f seconds', t(i)));
    grid on;
    
    pause(0.001); % Pause for animation effect
end;
% 
% % Bode plot
% figure;
% bode(sys);
% grid on;
% title('Bode Plot of the System');