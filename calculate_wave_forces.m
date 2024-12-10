function [F_total_x, F_total_y] = calculate_wave_forces(frequency, AmplitudeX, AmplitudeY, t)
    % Constants
    rho = 1025; % Seawater density (kg/m^3)
    g = 9.81; % Gravitational acceleration (m/s^2)
    V = 100; % Displaced volume of the floating body (m^3)
    Cd = 1.0; % Drag coefficient
    Ca = 0.5; % Added mass coefficient
    Cm = 1 + Ca; % Inertia coefficient

    % Pontoon dimensions
    length = 110; % Length of pier (m)
    width = 23; % Width of pier (m)
    height = 8.5; % Height of pier (m)

    V = length*width*height*5;

    % Cross-sectional area of the floating body (m^2)
    A = length * width;

    % Wave angular frequency (rad/s)
    omega = 2 * pi * frequency;

    % Water particle velocity (u)
    u = [AmplitudeX * cos(omega * t); AmplitudeY * sin(omega * t)]; % Assume 2D wave velocity field
    u_x = u(1, :); % x-direction velocity
    u_y = u(2, :); % y-direction velocity

    % Water particle acceleration (du/dt)
    du_x = AmplitudeX * -omega * sin(omega * t); % x-direction acceleration
    du_y = AmplitudeY * omega * cos(omega * t); % y-direction acceleration

    % Object velocity (v) is zero
    v_x = 0; % x-direction velocity
    v_y = 0; % y-direction velocity

    % Force calculations
    % (a) Froude-Krylov force
    Froude_Krylov_x = rho * V * du_x; % x-direction
    Froude_Krylov_y = rho * V * du_y; % y-direction

    % (b) Hydrodynamic mass force (Added mass)
    Hydrodynamic_mass_x = rho * Ca * V * (du_x - v_x); % x-direction
    Hydrodynamic_mass_y = rho * Ca * V * (du_y - v_y); % y-direction

    % (c) Drag force
    Drag_x = 0.5 * rho * Cd * A * (u_x - v_x) .* abs(u_x - v_x); % x-direction
    Drag_y = 0.5 * rho * Cd * A * (u_y - v_y) .* abs(u_y - v_y); % y-direction

    % Total dynamic force
    F_dynamic_x = Froude_Krylov_x + Hydrodynamic_mass_x + Drag_x;
    F_dynamic_y = Froude_Krylov_y + Hydrodynamic_mass_y + Drag_y;

    % Static buoyancy force
    F_static = rho * g * V % Acts vertically upwards

    % Total force (sum of dynamic and static components)
    F_total_x = F_dynamic_x; % No static component in x-direction
    F_total_y = F_dynamic_y + F_static; % Static buoyancy adds to y-direction force
    F_total_y = F_static; % Static buoyancy adds to y-direction force
end
% 
% rho = 1025; % Seawater density (kg/m^3)
% g = 9.81; % Gravitational acceleration (m/s^2)
% V = 100; % Displaced volume of the floating body (m^3)
% Cd = 1.0; % Drag coefficient
% Ca = 0.5; % Added mass coefficient
% Cm = 1 + Ca; % Inertia coefficient
% 
% % pier length, width, height
% length = 110;
% width = 23;
% height = 8.5;
% 
% waveFrequency = 0.2;
% t = linspace(0, 50, 500); % Time array (s)
% AmplitudeX = 5;
% AmplitudeY = 10;
% 
% A = length * width; % Cross-sectional area of the floating body (m^2)
% 
% omega = 2 * pi * waveFrequency; % Wave angular frequency (rad/s)
% 
% % Water particle velocity (u)
% u = [AmplitudeX * cos(omega * t); AmplitudeY * sin(omega * t)]; % Assume 2D wave velocity field
% u_x = u(1, :); % x-direction velocity
% u_y = u(2, :); % y-direction velocity
% 
% du_x = AmplitudeX * -omega * sin(omega * t); % x-direction acceleration
% du_y = AmplitudeY * omega * cos(omega * t); % y-direction acceleration
% 
% % Object velocity (v) is zero
% v_x = 0; % x-direction velocity
% v_y = 0; % y-direction velocity
% 
% % Force calculations
% % (a) Froude-Krylov force
% Froude_Krylov_x = rho * V * du_x; % x-direction
% Froude_Krylov_y = rho * V * du_y; % y-direction
% 
% % (b) Hydrodynamic mass force (Added mass)
% Hydrodynamic_mass_x = rho * Ca * V * (du_x - 0); % x-direction
% Hydrodynamic_mass_y = rho * Ca * V * (du_y - 0); % y-direction
% 
% % (c) Drag force
% Drag_x = 0.5 * rho * Cd * A * (u_x - v_x) .* abs(u_x - v_x); % x-direction
% Drag_y = 0.5 * rho * Cd * A * (u_y - v_y) .* abs(u_y - v_y); % y-direction
% 
% % Total dynamic force
% F_dynamic_x = Froude_Krylov_x + Hydrodynamic_mass_x + Drag_x;
% F_dynamic_y = Froude_Krylov_y + Hydrodynamic_mass_y + Drag_y;
% 
% % Static buoyancy force
% F_static = rho * g * V; % Acts vertically upwards
% 
% % Total force (sum of dynamic and static components)
% F_total_x = F_dynamic_x; % No static component in x-direction
% F_total_y = F_dynamic_y + F_static; % Static buoyancy adds to y-direction force
% 
% % Plot results
% figure;
% subplot(3, 1, 1);
% plot(t, F_dynamic_x, 'r', t, F_dynamic_y, 'b');
% xlabel('Time (s)');
% ylabel('Dynamic Force (N)');
% legend('F_{dynamic,x}', 'F_{dynamic,y}');
% title('Dynamic Forces');
% grid on;
% 
% subplot(3, 1, 2);
% plot(t, Froude_Krylov_x, 'g', t, Hydrodynamic_mass_x, 'm', t, Drag_x, 'k');
% xlabel('Time (s)');
% ylabel('Force Components (N)');
% legend('F_{Froude-Krylov}', 'F_{Hydrodynamic Mass}', 'F_{Drag}');
% title('Force Components in x-Direction');
% grid on;
% 
% subplot(3, 1, 3);
% plot(t, F_total_y, 'b');
% xlabel('Time (s)');
% ylabel('Total Force in y-Direction (N)');
% title('Total Force in Vertical Direction (Including Buoyancy)');
% grid on;