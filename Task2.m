clear; clc; close all;
%% Preparatory settings
res_fld = 'results';
NACA_4415 = load(fullfile(res_fld, 'XFOIL_NACA_4415.mat')).NACA_4415;

%% Simulation settings
alpha = [0 5 10];
N=100;

%% Lifting Line calculations
AR = 4:2:10;
AR = [AR, 1e+5];  % Add 'infinite' aspect ratio
b = 1;
for i = 1:numel(AR)
    wings_rect(i) = struct('AR', AR(i), 'b', b, 'wing', ...
                           RectangularWing(AR(i), b));
end

for i = 1:numel(AR)
    % Calculate coefficients of lifting line theory
    [y, theta] = wings_rect(i).wing.generate_coordinates(N);
    A = LiftingLine.solve_coeffs(wings_rect(i).wing, y, theta, alpha, N, ...
        NACA_4415.m_0, NACA_4415.alpha_L0);
    % Calculate spanwise parameters
    [alpha_i, C_l, C_di, Gamma] = ...
        LiftingLine.calc_lift_drag_sections(wings_rect(i).wing, y, theta, A);
    
    alpha_g = alpha - alpha_i;  % Geometric angle of attack
    C_d_fric = interp1(NACA_4415.xfoil_res.alpha, ...
                       NACA_4415.xfoil_res.C_d, ...
                       alpha_g);
    wings_rect(i).LL_res = struct('y', y, 'theta', theta, 'A', A, ...
                                 'alpha', alpha,...
                                 'alpha_i', alpha_i, 'C_l', C_l, ...
                                 'C_di', C_di, 'C_d_fric', C_d_fric, ...
                                 'Gamma', Gamma);
end

save(fullfile(res_fld, 'T2_wings_rect.mat'), 'wings_rect');

%% Plots
% Selection
savefigs = true;
plot_Gamma = true;
plot_alpha = true;
plot_C_l = false;
plot_C_di = false;

% Settings
cols = ["#0072BD", "#D95319", "#EDB120", "#77AC30", "#80B3FF"];  % Colors of the lines
markers = ["+", "*", "o", "diamond", "none"];  % Markers for the four methods
ms = [4.5, 4.5, 4.5, 4.5, 4.5];  % Marker size for the plots of the four methods
lw = [1, 1, 1, 1, 1.5];  % Linewidth for the lines of the four methods
ax_col = [0.2, 0.2, 0.2];  % Color of accented axes
ax_lw = 1.5;  % Line width of accented axes
fs = 16;  % Plot font size
fig_count = 0;

%% Preparation
% Create export directory if it doesn't exist
exp_fld = 'plots';
if ~exist(exp_fld, 'dir')
    mkdir(exp_fld);
end

%Plot Gamma vs y
if plot_Gamma
    for i_a = 1:numel(alpha)
        % Create plot
        figure(i_a+fig_count);
        cla; hold on; grid on;
        colororder(cols);
        ax = gca;
            
        % Highlight y=0 grid line
        y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                     HandleVisibility='off'); % Thick vertical line at x=0
    
        % Plot Gamma curves 
        for i = 1:numel(AR)
            if AR(i)<10000
                disp_name = sprintf('$AR=%d$', AR(i));
            else
                disp_name = '$AR=\infty$';
            end
            plot(wings_rect(i).LL_res.y/wings_rect(i).b, ...
                 wings_rect(i).LL_res.Gamma(:,i_a), ...
                 LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
                 DisplayName=disp_name);
        end
        hold off; 

        % Configure limits and ticks
        ylim('auto');
        xticks(-.5:.25:.5);
        xlim(ax, [-.5, .5]);
    
        % Plot labels
        set(gcf,'Color','White');
        set(ax,'FontSize',fs);
        legend('Location', 'south', 'Interpreter', 'latex')
        xlabel('$y/b$', 'Interpreter', 'latex');
        ylabel('$\Gamma$', 'Interpreter', 'latex');
        set(ax, 'TickLabelInterpreter', 'latex');
    end
else
    disp('Gamma vs y not plotted')
end

fig_count = fig_count + numel(alpha);

%Plot alpha_i vs y
if plot_alpha
    for i_a = 1:numel(alpha)
        % Create plot
        figure(i_a+fig_count);
        cla; hold on; grid on;
        colororder(cols);
        ax = gca;
            
        % Highlight y=0 grid line
        y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                     HandleVisibility='off'); % Thick vertical line at x=0
    
        % Plot alpha curves 
        plt = [];  % Ensure plot list is empty
        for i = 1:numel(AR)
            if AR(i)<10000
                disp_name = sprintf('$AR=%d$', AR(i));
            else
                disp_name = '$AR=\infty$';
            end
            plt(i) = plot(wings_rect(i).LL_res.y/wings_rect(i).b, ...
                          wings_rect(i).LL_res.alpha_i(:,i_a), ...
                          LineWidth=lw(i), ...
                          Marker=markers(i), MarkerSize=ms(i), ...
                          DisplayName=disp_name);
        end
        uistack(plt(end), 'bottom');
        hold off; 
    
        % % Configure limits and ticks
        ylim('auto');
        xticks(-.5:.25:.5);
        xlim(ax, [-.5, .5]);
    
        % Plot labels
        set(gcf,'Color','White');
        set(ax,'FontSize',fs);
        legend(plt, 'Location', 'north', 'Interpreter', 'latex')
        xlabel('$y/b$', 'Interpreter', 'latex');
        ylabel('$\alpha_i\:[^{\circ}]$', 'Interpreter', 'latex');
        set(ax, 'TickLabelInterpreter', 'latex');
    
        % Save figure
            if savefigs
                exp_name = fullfile(exp_fld, ...
                    sprintf('alpha_i_vs_y_AoA=%d.pdf', alpha(i_a)));
                exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                    'BackgroundColor', 'none', 'Resolution', 300);
            end
    end
else
    disp('C_l vs alpha not plotted')
end

fig_count = fig_count + numel(alpha);