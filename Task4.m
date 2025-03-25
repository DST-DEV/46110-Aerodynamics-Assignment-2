clear; clc; close all;
%% Preparatory settings
res_fld = 'results';
NACA_4415 = load(fullfile(res_fld, 'XFOIL_NACA_4415.mat')).NACA_4415;

%% Simulation settings
alpha_g0 = 4;  % Midpoint geometric angle of attack
alpha_gt = 0:2:8;  % Tip geometric angle of attack
U_0 = 1;
N=100;

%% Lifting Line calculations
AR = 6;
b = 1;
wing = RectangularWing(AR, b);
wings_rect = struct('AR', AR, 'b', b, 'alpha_g0', alpha_g0, ...
                       'alpha_gt', alpha_gt, 'wing', wing);

% Calculate coordinates and spanwise geometric angle of attack
% distribution
[y, theta] = wing.generate_coordinates(N);
alpha_g = alpha_g0 + (alpha_gt-alpha_g0).*abs(y'/b*2);

% Calculate coefficients of lifting line theory
A = LiftingLine.solve_coeffs(wing, y, theta, alpha_g, ...
    NACA_4415.m_0, NACA_4415.alpha_L0);

% Calculate spanwise parameters
[alpha_i, C_l, C_di, Gamma] = ...
    LiftingLine.calc_lift_drag_sections(wing, y, theta, A);

%Non-dimensionalize Gamma distribution
Gamma_nd = Gamma./(U_0*wing.chord_length(y)'./wing.c_root);

wings_rect.LL_res = struct('y', y, 'theta', theta, 'A', A, ...
                             'alpha_g', alpha_g,...
                             'alpha_i', alpha_i, 'C_l', C_l, ...
                             'C_di', C_di,  ...
                             'Gamma', Gamma, 'Gamma_nd', Gamma_nd);

save(fullfile(res_fld, 'T4_wings_rect_twist.mat'), 'wings_rect');

%% Plots
% Selection
savefigs = true;
plot_Gamma = true;
plot_alpha = true;
plot_C_l = true;
plot_C_di = true;

% Settings
cols = ["#0072BD", "#D95319", "#EDB120", "#77AC30", "#80B3FF"];  % Colors of the lines
markers = ["+", "*", "o", "diamond", "v"];  % Markers for the four methods
ms = [4.5, 4.5, 4.5, 4.5, 4.5];  % Marker size for the plots of the four methods
lw = [1, 1, 1, 1, 1];  % Linewidth for the lines of the four methods
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
    % Create plot
    figure(fig_count + 1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot Gamma curves 
    for i = 1:numel(alpha_gt)
        plot(wings_rect.LL_res.y/wings_rect.b*2, ...
             wings_rect.LL_res.Gamma_nd(:,i), ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=sprintf('$\\alpha_{g,tip}=%d$', alpha_gt(i)));
    end
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'south', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$\Gamma$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T4_Gamma_vs_y.pdf');
        exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('Gamma vs y not plotted')
end

fig_count = fig_count + 1;

%Plot alpha_i vs y
if plot_alpha
    % Create plot
    figure(fig_count + 1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot Gamma curves 
    for i = 1:numel(alpha_gt)
        plot(wings_rect.LL_res.y/wings_rect.b*2, ...
             wings_rect.LL_res.alpha_i(:,i), ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=[sprintf('$\\alpha_{g,tip}=%d$', alpha_gt(i))]);
    end
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'north', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$\alpha_i\:[^{\circ}]$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T4_alpha_i_vs_y.pdf');
        exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('induced alpha vs y not plotted')
end

fig_count = fig_count + 1;

%Plot C_l vs y
if plot_C_l
    % Create plot
    figure(fig_count + 1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot C_l curves 
    for i = 1:numel(alpha_gt)
        plot(wings_rect.LL_res.y/wings_rect.b*2, ...
             wings_rect.LL_res.C_l(:,i), ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=sprintf('$\\alpha_{g,tip}=%d$', alpha_gt(i)));
    end
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'south', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$C_l$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T4_C_l_vs_y.pdf');
        exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('C_l vs y not plotted')
end

fig_count = fig_count + 1;

%Plot C_d vs y
if plot_C_di
    % Create plot
    figure(fig_count + 1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot C_d curves 
    for i = 1:numel(alpha_gt)
        plot(wings_rect.LL_res.y/wings_rect.b*2, ...
             wings_rect.LL_res.C_di(:,i), ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=sprintf('$\\alpha_{g,tip}=%d$', alpha_gt(i)));
    end
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'north', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$C_{d_i}$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T4_C_di_vs_y.pdf');
        exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('C_d_induced vs y not plotted')
end

fig_count = fig_count + 1;