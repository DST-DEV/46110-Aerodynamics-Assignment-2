clear; clc; close all;
%% Preparatory settings
res_fld = 'results';
NACA_4415 = load(fullfile(res_fld, 'XFOIL_NACA_4415.mat')).NACA_4415;

%% Simulation settings
alpha = 6;
U_0 = 1;
N=201;

%% Lifting Line calculations
AR = 6;
TR = .2:.2:1;
b = 1;
for i = 1:numel(TR)
    wings_taper(i) = struct('AR', AR, 'TR', TR(i), 'b', b, ...
                            'wing', TaperedWing(AR, TR(i), b));
end
for i = 1:numel(TR)
    wing_i = wings_taper(i).wing;
    % Calculate spanwise coordinates 
    [y, theta] = wing_i.generate_coordinates(N);

    % Calculate coefficients of lifting line theory
    A = LiftingLine.solve_coeffs(wing_i, y, theta, alpha, ...
        NACA_4415.m_0, NACA_4415.alpha_L0);

    % Calculate spanwise parameters
    [alpha_i, C_l, C_di, Gamma] = ...
        LiftingLine.calc_lift_drag_sections(wing_i, y, theta, A);

    %Non-dimensionalize Gamma distribution
    Gamma_nd = Gamma./(U_0*wing_i.c_mean);
    
    wings_taper(i).LL_res = struct('y', y, 'theta', theta, 'A', A, ...
                                 'alpha', alpha,...
                                 'alpha_i', alpha_i, 'C_l', C_l, ...
                                 'C_di', C_di, 'Gamma', Gamma, ...
                                 'Gamma_nd', Gamma_nd);
end

save(fullfile(res_fld, 'T3_wings_tapered.mat'), 'wings_taper');

%% Plots
% Selection
savefigs = true;
plot_Gamma = true;
plot_alpha = true;
plot_C_l = true;
plot_C_di = true;

% Settings
cols = ["#0072BD", "#D95319", "#EDB120", "#77AC30", "#80B3FF"];  % Colors of the lines
markers = ["none", "none", "none", "none", "none"];  % Markers for the four methods
ms = [4.5, 4.5, 4.5, 4.5, 4.5];  % Marker size for the plots of the four methods
lw = [1.5, 1.5, 1.5, 1.5, 1.5];  % Linewidth for the lines of the four methods
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
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot Gamma curves 
    for i = 1:numel(TR)
        plot(wings_taper(i).LL_res.y/wings_taper(i).b*2, ...
             wings_taper(i).LL_res.Gamma_nd, ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=sprintf('$TR=%.1f$', TR(i)));
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
    ylabel('$\Gamma/(U_\infty\bar{c})$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T3_Gamma_vs_y.pdf');
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
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot alpha curves 
    plt = [];  % Ensure plot list is empty
    for i = 1:numel(TR)
        plt(i) = plot(wings_taper(i).LL_res.y/wings_taper(i).b*2, ...
                      wings_taper(i).LL_res.alpha_i, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=sprintf('$TR=%.1f$', TR(i)));
    end
    hold off; 

    % % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend(plt, 'Location', 'north', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$\alpha_i\:[^{\circ}]$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T3_alpha_i_vs_y.pdf');
    exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
        'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('Induced alpha vs y not plotted')
end

fig_count = fig_count + 1;

%Plot C_l vs y
if plot_C_l
    % Create plot
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot alpha curves 
    plt = [];  % Ensure plot list is empty
    for i = 1:numel(TR)
        plt(i) = plot(wings_taper(i).LL_res.y/wings_taper(i).b*2, ...
                      wings_taper(i).LL_res.C_l, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=sprintf('$TR=%.1f$', TR(i)));
    end
    hold off; 

    % % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend(plt, 'Location', 'south', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$C_l$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T3_C_l_vs_y.pdf');
    exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
        'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('C_l vs y not plotted')
end

fig_count = fig_count + 1;

%Plot C_di vs y
if plot_C_di
    % Create plot
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot alpha curves 
    plt = [];  % Ensure plot list is empty
    for i = 1:numel(TR)
        plt(i) = plot(wings_taper(i).LL_res.y/wings_taper(i).b*2, ...
                      wings_taper(i).LL_res.C_di, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=sprintf('$TR=%.1f$', TR(i)));
    end
    hold off; 

    % % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend(plt, 'Location', 'south', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$C_{d_i}$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'T3_C_l_vs_y.pdf');
    exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
        'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('C_di vs y not plotted')
end

fig_count = fig_count + 1;