clear; clc; close all;
%% Preparatory settings
res_fld = 'results';
NACA_4415 = load(fullfile(res_fld, 'XFOIL_NACA_4415.mat')).NACA_4415;

%% Simulation settings
alpha = -6:.25:10;
U_0 = 1;
N=101;

%% Lifting Line calculations
AR = 4:2:10;
AR = [AR, 1e+50];  % Add 'infinite' aspect ratio
b = 1;
for i = 1:numel(AR)
    wings_ell(i) = struct('AR', AR(i), 'b', b, 'wing', ...
                          EllipticWing(AR(i), b));
end

for i = 1:numel(AR)
    % Calculate spanwise coordinates 
    [y, theta] = wings_ell(i).wing.generate_coordinates(N);

    % Calculate coefficients of lifting line theory
    A = LiftingLine.solve_coeffs(wings_ell(i).wing, y, theta, alpha, ...
        NACA_4415.m_0, NACA_4415.alpha_L0);

    % Calculate global lift and induced drag coefficients
    [C_l_tot, C_di_tot] = LiftingLine.calc_lift_drag_wing(wings_ell(i).wing, A);

    % Calculate spanwise parameters
    % Note: alpha_i should be constant along the wing span
    [alpha_i, C_l, C_di, Gamma] = ...
        LiftingLine.calc_lift_drag_sections(wings_ell(i).wing, y, theta, A);
    
    % Calculate friction drag
    alpha_g = alpha - alpha_i;  % Geometric angle of attack
    C_d_fric = interp1(NACA_4415.xfoil_res.alpha, ...
                       NACA_4415.xfoil_res.C_d, ...
                       alpha_g);
    c = wings_ell(i).wing.chord_length(y);
    D = -trapz(y, c'.*C_d_fric, 1);
    C_d_fric_tot = D/wings_ell(i).wing.S;

    %Non-dimensionalize Gamma distribution
    Gamma_nd = Gamma./(U_0*wings_ell(i).wing.c_mean);

    wings_ell(i).LL_res = struct('y', y, 'theta', theta, 'A', A, 'c', c, ...
                                 'alpha', alpha,...
                                 'C_l_tot', C_l_tot, ...
                                 'C_di_tot', C_di_tot, ...
                                 'alpha_i', alpha_i, 'C_l', C_l, ...
                                 'C_di', C_di, 'C_d_fric', C_d_fric, ...
                                 'D', D, 'C_d_fric_tot', C_d_fric_tot,...
                                 'Gamma', Gamma, 'Gamma_nd', Gamma_nd);
end

save(fullfile(res_fld, 'T1_wings_ell.mat'), 'wings_ell');

%% Plots
% Selection
savefigs = true;
plot_Gamma = true; gamma_aoa = 5;
plot_alpha_i = true;
plot_C_l = true;
plot_C_di = true;
plot_C_D = true;

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
    ind = find(alpha == gamma_aoa);

    % Create plot
    figure(fig_count+1);
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
        plot(wings_ell(i).LL_res.y/wings_ell(i).b*2, ...
             wings_ell(i).LL_res.Gamma_nd(:,ind), ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=disp_name);
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
else
    disp('Gamma vs y not plotted')
end

fig_count = fig_count + 1;

%Plot alpha_i vs alpha
if plot_alpha_i
    % Create plot
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot C_l curves 
    for i = 1:numel(AR)
        if AR(i)<10000
            disp_name = sprintf('$AR=%d$', AR(i));
        else
            disp_name = '$AR=\infty$';
        end
        plt(i) = plot(wings_ell(i).LL_res.alpha, ...
                      wings_ell(i).LL_res.alpha_i(1,:), ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=disp_name);
    end
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(min(alpha):2:max(alpha));
    xlim(ax, [min(alpha), max(alpha)]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'northwest', 'Interpreter', 'latex')
    xlabel('AoA $[^{\circ}]$', 'Interpreter', 'latex');
    ylabel('$\alpha_i\:[^{\circ}]$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, 'T1_alpha_i_vs_alpha.pdf');
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
else
    disp('alpha_i vs alpha not plotted')
end

fig_count = fig_count + 1;

%Plot C_L vs alpha
if plot_C_l
    % Create plot
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight x=0 and y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0
    x_ax = yline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at y=0

    % Plot C_L curves 
    for i = 1:numel(AR)
        if AR(i)<10000
            disp_name = sprintf('$AR=%d$', AR(i));
        else
            disp_name = '$AR=\infty$';
        end
        plt(i) = plot(wings_ell(i).LL_res.alpha, ...
                      wings_ell(i).LL_res.C_l_tot, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=disp_name);
    end
    plt(end+1) = plot(alpha, deg2rad(alpha-NACA_4415.alpha_L0)*NACA_4415.m_0, ...
                  LineWidth=1, color="k", LineStyle = "-.",...
                  DisplayName='Linear Lift Slope');
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(min(alpha):2:max(alpha));
    xlim(ax, [min(alpha), max(alpha)]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'northwest', 'Interpreter', 'latex')
    xlabel('AoA $[^{\circ}]$', 'Interpreter', 'latex');
    ylabel('$C_L$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, 'T1_C_l_vs_alpha.pdf');
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
else
    disp('C_l vs alpha not plotted')
end

fig_count = fig_count + 1;

%Plot C_Di vs alpha
if plot_C_di
    % Create plot
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
    
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot C_Di curves 
    plt = [];  % Ensure plot list is empty
    for i = 1:numel(AR)
        if AR(i)<10000
            disp_name = sprintf('$AR=%d$', AR(i));
        else
            disp_name = '$AR=\infty$';
        end
        plt(i) = plot(wings_ell(i).LL_res.alpha, ...
                      wings_ell(i).LL_res.C_di_tot, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=disp_name);
    end
    uistack(plt(end), 'bottom');
    hold off; 

    % Configure limits and ticks
    ylim(ax, [-.005, .1]);
    xticks(min(alpha):2:max(alpha));
    xlim(ax, [min(alpha), max(alpha)]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend(plt, 'Location', 'northwest', 'Interpreter', 'latex')
    xlabel('AoA $[^{\circ}]$', 'Interpreter', 'latex');
    ylabel('$C_{D_i}$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, 'T1_C_di_vs_alpha.pdf');
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
else
    disp('C_d induced vs alpha not plotted')
end

fig_count = fig_count + 1;

%Plot C_D vs alpha
if plot_C_D
    % Create plot
    figure(fig_count+1);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
    
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot C_d curves 
    plt = [];  % Ensure plot list is empty
    for i = 1:numel(AR)
        if AR(i)<10000
            disp_name = sprintf('$AR=%d$', AR(i));
        else
            disp_name = '$AR=\infty$';
        end
        plt(i) = plot(wings_ell(i).LL_res.alpha, ...
                      wings_ell(i).LL_res.C_d_fric_tot, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=disp_name);
    end
    plt(end+1) = plot(NACA_4415.xfoil_res.alpha, ...
                      NACA_4415.xfoil_res.C_d, ...
                      LineWidth=1, color="k", LineStyle = "-.", ...
                      DisplayName='NACA 4415');
    hold off; 

    % Configure limits and ticks
    xticks(min(alpha):2:max(alpha));
    xlim(ax, [min(alpha), max(alpha)]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend(plt, 'Location', 'northwest', 'Interpreter', 'latex')
    xlabel('AoA $[^{\circ}]$', 'Interpreter', 'latex');
    ylabel('$C_{D_{fric}}$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, 'T1_C_D_vs_alpha.pdf');
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
else
    disp('C_D vs alpha not plotted')
end

fig_count = fig_count + 1;