%
% MATLAB Code for 1D Convergent-Divergent Nozzle Analysis
%

% Header
close all
clearvars
clc

% Load geometry data
load("nozzle_geometry.mat")
geom = geometry_data;

% Input Parameters
gamma = 1.4;       % Specific heat ratio (-)
R = 287.1;         % Specific gas constant (J/ kg / K)
p0 = 200000;       % Stagnation inlet pressure (Pa)
T0 = 300;          % Stagnation inlet temperature (K)

% Geometrical specifications
x_throat = 50 * 1e-3;                           % Throat axial position (m)
L_throat = 5 * 1e-3;                            % length of throat (m)
L_div = 120 * 1e-3 - (x_throat + L_throat);     % outlet axial position (m)
L_nozzle = x_throat + L_throat + L_div;         % Nozzle total length (m)
thickness = 1 * 1e-3;                           % Thickness (m)
x = linspace(0, L_nozzle, 3e3);                 % Axial positions (m)
x_geom = geom.x;                                % vector of known axial positions (m)
h_geom = geom.h;                                % vector of known nozzle wall positions (m)

% Calculation of Area at desired axial positions
A = interp1(geom.x,geom.h,x) * thickness;
A_throat = interp1(geom.x,geom.h,x_throat) * thickness;
A_out = interp1(geom.x,geom.h,L_nozzle) * thickness;

% alternative geometry distributions
% A = area_distribution(x,A_in,A_out,A_throat,x_throat,L_nozzle); % Area distribution (two parabolic profiles)
% A = area_distribution(x,D_in,D_out,D_throat,x_throat,L_nozzle); % Area distribution (two third-degree profiles)

% Calculation of back pressure that makes throat be chocked
M_throat = 1;   % (-) - Suppose the throat is choked
f0= @(M_out) A_throat/A_out - M_out/M_throat * ...
    (1 + (gamma - 1) / 2 * M_throat^2)^((gamma + 1) / (2 * (gamma - 1))) /...
    (1 + (gamma - 1) / 2 * M_out^2)^((gamma + 1) / (2 * (gamma - 1)));  % relationship between area distribution and Mach
M_out = fzero(f0, [0.001 1-1e-8]);  % (-) numerical resolution to calculate Mach
pb_max = p0 * (1 + (gamma-1)/2 * M_out^2)^(-gamma/(gamma-1)); % (Pa) maximum back pressure that chokes the throat

% Back pressure
pb = 1.5 * 1e5; % Set back pressure (Pa)
% pb = pout_sub; % Set back pressure (Pa)
if pb >= pb_max % to avoid simulating unsupported configurations (i.e, throat not chocked) 
    error('The code does not support back pressure values that DO NOT choke the throat. Use pb lower than pb_max.')
end

% Initialize variables (for speed)
M = nan(size(x));     % Mach number
p = nan(size(x));     % Pressure
T = nan(size(x));     % Temperature
rho = nan(size(x));   % Density
shock_detected = false; % Flag for shock detection

% Calculate flow properties
for i = 1:length(x)
    if abs(A(i) - A_throat)/A_throat < 1e-8  % checks whether we are in the throat
        % Throat condition: M = 1
        M(i) = 1;
    else
        if x(i) <= x_throat
            % Subsonic flow in convergent section
            % Area-Mach relation (solved numerically)
            func = @(M) A(i)/A_throat - 1/M * (2/(gamma+1) * (1 + (gamma-1) / 2 * M^2))^((gamma+1) / (2*(gamma-1)));
            M(i) = fzero(func, [0.001 1-1e-8]);
        elseif x(i) > x_throat && shock_detected == false && pb >= pb_max
            % Subsonic flow in divergent section before shock
            % Area-Mach relation (solved numerically)
            func = @(M) A(i)/A_throat - 1/M * (2/(gamma+1) * (1 + (gamma-1) / 2 * M^2))^((gamma+1) / (2*(gamma-1)));
            M(i) = fzero(func, [0.001 1-1e-8]);
        elseif x(i) > x_throat && shock_detected == false
            % Supersonic flow in divergent section before shock
            % Area-Mach relation (solved numerically)
            func = @(M) A(i)/A_throat - 1/M * (2/(gamma+1) * (1 + (gamma-1) / 2 * M^2))^((gamma+1) / (2*(gamma-1)));
            M(i) = fzero(func, [1+1e-8 5]);
        elseif x(i) > x_throat && shock_detected == true
            % Subsonic flow in divergent section after shock
            % Area-Mach relation (solved numerically)
            func_after_shock= @(M) A(i)/A_shock - M_after_shock/M * ...
                (1 + (gamma - 1) / 2 * M^2)^((gamma + 1) / (2*(gamma - 1))) /...
            (1 + (gamma - 1) / 2 * M_after_shock^2)^((gamma + 1) / (2*(gamma - 1)));
            M(i) = fzero(func_after_shock, [0.001 1-1e-8]);
        else
            error('This should not happen')
        end
    end

    % Isentropic relations
    if shock_detected == false
        % Before shock p0 and T0 are used
        p(i) = p0 * (1 + (gamma-1)/2 * M(i)^2)^(-gamma/(gamma-1));
        T(i) = T0 * (1 + (gamma-1)/2 * M(i)^2)^-1;
    else
        % After shock p2 and T2 are used
        p(i) = p20 * (1 + (gamma-1)/2 * M(i)^2)^(-gamma/(gamma-1));
        T(i) = T20 * (1 + (gamma-1)/2 * M(i)^2)^-1;
    end
    % Density calculation
    rho(i) = p(i) / (R * T(i));

    % Identify shock location and recalculate total pressure and
    % temperature
    if x(i) >= x_throat + L_throat && p(i) < pb && shock_detected == false && pb < pb_max
        
        % we suppose to have the shock at the i-th x value
        A_shock_guess = A(i);
        
        % Normal shock relations at shock location
        M1_guess = M(i); % Upstream Mach number
        M2_guess = sqrt( ...
            (2 + (gamma-1) * M1_guess^2) / ...
            (2 * gamma * M1_guess^2 - (gamma-1))); % Downstream Mach number
        
        % Recalculate total pressure and temperature after the shock
        p20_guess = p0 * (...
        ((gamma + 1) * M1_guess^2) / ((gamma - 1) * M1_guess^2 + 2))^(gamma / (gamma-1)) * ...
        ((gamma + 1) / (2 * gamma * M1_guess^2 - (gamma-1)))^(1 / (gamma-1)); % Downstream total pressure
        T20_guess = T0;
        
        % Find Mach at the exit 
        func_out= @(M_out) A_out/A_shock_guess - M2_guess/M_out * ...
            (1 + (gamma - 1) / 2 * M_out^2)^((gamma + 1) / (2*(gamma - 1))) /...
            (1 + (gamma - 1) / 2 * M2_guess^2)^((gamma + 1) / (2*(gamma - 1)));
        M_out_guess = fzero(func_out, [0.001 1-1e-8]); % Subsonic flow after the shock
        
        % Calculate the p_out in case the shock has occurred
        p_out_guess = p20_guess * (1 + (gamma-1)/2 * M_out_guess^2)^(-gamma/(gamma-1));

        % Determine if the shock actually occurred
        if p_out_guess >= pb     
           % no shock - nothing happens
        else                     
            % the shock actually occured
            shock_detected = true;
            p20 = p20_guess;
            T20 = T20_guess;
            A_shock = A_shock_guess;
            M_after_shock = M2_guess;
        end
    end
end

% Plot results
figure
tiledlayout(2,2,"TileSpacing","compact","Padding","compact")

nexttile
plot(x, A, "LineWidth",2)
xlabel('x (m)')
ylabel('Area (m2)')
title('Area vs. Position')
set(gca,"FontSize",16,"FontName","Calibri")
ylim([0 inf])
grid on

nexttile
plot(x, M, "LineWidth",2)
xlabel('x (m)')
ylabel('Mach Number')
title('Mach Number vs. Position')
set(gca,"FontSize",16,"FontName","Calibri")
grid on

nexttile
plot(x, p, "LineWidth",2)
yline(pb,"LineWidth",2)
xlabel('x (m)')
ylabel('Pressure (Pa)')
title('Pressure vs. Position')
set(gca,"FontSize",16,"FontName","Calibri")
grid on
legend("Pressure","Back Pressure",'location','southwest')

nexttile
plot(x, T, "LineWidth",2)
xlabel('x (m)')
ylabel('Temperature (K)')
title('Temperature vs. Position')
set(gca,"FontSize",16,"FontName","Calibri")
grid on

% Comparison with CFD data
% CFD without friction and viscosity
load("cfd_data_freeslip.mat")   % load CFD data

figure
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")

nexttile
idx_toplot = pressure_data.x_m <= 120 * 1e-3;
plot(pressure_data.x_m(idx_toplot),pressure_data.P110000_Pa(idx_toplot) / 1e5,"LineWidth",2)
hold on
plot(pressure_data.x_m(idx_toplot),pressure_data.P135000_Pa(idx_toplot) / 1e5,"LineWidth",2)
plot(pressure_data.x_m(idx_toplot),pressure_data.P150000_Pa(idx_toplot) / 1e5,"LineWidth",2)
plot(x,p / 1e5,'LineWidth',2,"LineWidth",2)
xlabel('x (m)')
ylabel('Pressure (bar)')
grid on
legend("CFD (p_{b} = 1.10 bar)","CFD (p_{b} = 1.35 bar)","CFD (p_{b} = 1.50 bar)",...
    string(['Model (p_{b} = ',num2str(pb / 1e5,3),' bar)']),...
    'Location','southwest')
set(gca,"FontSize",18,"FontName","Calibri")

nexttile
idx_toplot = mach_data.x_m <= 120 * 1e-3;
plot(mach_data.x_m(idx_toplot),mach_data.P110000_Pa(idx_toplot),"LineWidth",2)
hold on
plot(mach_data.x_m(idx_toplot),mach_data.P135000_Pa(idx_toplot),"LineWidth",2)
plot(mach_data.x_m(idx_toplot),mach_data.P150000_Pa(idx_toplot),"LineWidth",2)
plot(x,M,'LineWidth',2)
xlabel('x (m)')
ylabel('Mach (-)')
grid on
legend("CFD (p_{b} = 1.10 bar)","CFD (p_{b} = 1.35 bar)","CFD (p_{b} = 1.50 bar)",...
    string(['Model (p_{b} = ',num2str(pb / 1e5,3),' bar)']),...
    'Location','northwest')
set(gca,"FontSize",18,"FontName","Calibri")

% CFD with friction and viscosity
load("cfd_data_noslip.mat")   % load CFD data

figure
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")

nexttile
idx_toplot = pressure_data.x_m <= 120 * 1e-3;
plot(pressure_data.x_m(idx_toplot),pressure_data.P90000_Pa(idx_toplot) / 1e5,"LineWidth",2)
hold on
plot(pressure_data.x_m(idx_toplot),pressure_data.P110000_Pa(idx_toplot) / 1e5,"LineWidth",2)
plot(pressure_data.x_m(idx_toplot),pressure_data.P185000_Pa(idx_toplot) / 1e5,"LineWidth",2)
plot(x,p / 1e5,'LineWidth',2,"LineWidth",2)
xlabel('x (m)')
ylabel('Pressure (bar)')
grid on
legend("CFD (p_{b} = 0.90 bar)","CFD (p_{b} = 1.10 bar)","CFD (p_{b} = 1.85 bar)",...
    string(['Model (p_{b} = ',num2str(pb / 1e5,3),' bar)']),...
    'Location','southwest')
set(gca,"FontSize",18,"FontName","Calibri")

nexttile
idx_toplot = mach_data.x_m <= 120 * 1e-3;
plot(mach_data.x_m(idx_toplot),mach_data.P90000_Pa(idx_toplot),"LineWidth",2)
hold on
plot(mach_data.x_m(idx_toplot),mach_data.P110000_Pa(idx_toplot),"LineWidth",2)
plot(mach_data.x_m(idx_toplot),mach_data.P185000_Pa(idx_toplot),"LineWidth",2)
plot(x,M,'LineWidth',2)
xlabel('x (m)')
ylabel('Mach (-)')
grid on
legend("CFD (p_{b} = 0.90 bar)","CFD (p_{b} = 1.10 bar)","CFD (p_{b} = 1.85 bar)",...
    string(['Model (p_{b} = ',num2str(pb / 1e5,3),' bar)']),...
    'Location','northwest')
set(gca,"FontSize",18,"FontName","Calibri")



% Alternative area distributions

% Area distribution
% function [A] = area_distribution(x,D_in,D_out,D_throat,x_throat,L_nozzle)
% %
% A = nan(size(x));   % area initialised for speed
% for i = 1 : length(x)
%     if x(i) <= x_throat
%         % convergent part
%         z = x(i) / x_throat;
%         D = D_in + (D_in - D_throat) * (-10 * z^3 + 15 * z^4 - 6 * z^5);
%         A(i) = pi/4 * D^2;
%     else
%         % divergent part
%         z = (x(i) - x_throat) / (L_nozzle - x_throat);
%         D = D_throat + (D_throat - D_out) * (-10 * z^3 + 15 * z^4 - 6 * z^5);
%         A(i) = pi/4 * D^2;
%     end
% end
%     %
% end

% % Area distribution
% function [A] = area_distribution(x,A_in,A_out,A_throat,x_throat,L_nozzle)
% %
% A = nan(size(x));   % area initialised for speed
% for i = 1 : length(x)
%     if x(i) <= x_throat
%         % convergent part
%         A(i) = (A_in - A_throat) / x_throat^2 * (x(i) - x_throat)^2 + A_throat;
%     else
%         % divergent part
%         A(i) = (A_out - A_throat) / (L_nozzle - x_throat)^2 * (x(i) - x_throat).^2 + A_throat;
%     end
% end
%     %
% end
%
%