clear; clc;
%% Preparatory settings
res_fld = 'results';
NACA_4415 = load(fullfile(res_fld, 'XFOIL_NACA_4415.mat')).NACA_4415;

%% Simulation settings
alpha = -4:.1:10;
i_aa = find(alpha == 0 | alpha == 5 | alpha == 10); % Find indices for sub-taks a
N=101;

%% Lifting Line calculations
AR = 4:2:10;
AR = [AR, 1e+5];  % Add 'infinite' aspect ratio
b = 1;
for i = 1:numel(AR)
    wings_rect(i) = struct('AR', AR(i), 'b', b, 'wing', ...
                           RectangularWing(AR(i), b));
end

for i = 1:numel(AR)
    % Calculate spanwise coordinates 
    [y, theta] = wings_rect(i).wing.generate_coordinates(N, 0);
    figure(9)
    plot (y)
    figure(10)
    plot (theta)

    % Calculate coefficients of lifting line theory
    A = LiftingLine.solve_coeffs(wings_rect(i).wing, y, theta, alpha, ...
        NACA_4415.m_0, NACA_4415.alpha_L0);
    % Calculate global lift and induced drag coefficients
    [C_l_tot, C_di_tot] = LiftingLine.calc_lift_drag_wing(wings_rect(i).wing, A);

    % Calculate spanwise parameters
    [alpha_i, C_l, C_di, Gamma] = ...
        LiftingLine.calc_lift_drag_sections(wings_rect(i).wing, y, theta, A);
    
    alpha_eff = alpha - alpha_i;  % Effective angle of attack
    C_d_fric = interp1(NACA_4415.xfoil_res.alpha, ...
                       NACA_4415.xfoil_res.C_d, ...
                       alpha_eff);
    wings_rect(i).LL_res = struct('y', y, 'theta', theta, 'A', A, ...
                                 'alpha', alpha,...
                                 'C_l_tot', C_l_tot, ...
                                 'C_di_tot', C_di_tot, ...
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
    for i_a = 1:length(i_aa)
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
            plot(wings_rect(i).LL_res.y/wings_rect(i).b*2, ...
                 wings_rect(i).LL_res.Gamma(:,i_aa(i_a)), ...
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
        ylabel('$\Gamma$', 'Interpreter', 'latex');
        set(ax, 'TickLabelInterpreter', 'latex');
    end
else
    disp('Gamma vs y not plotted')
end

fig_count = fig_count + length(i_aa);

%Plot alpha_i vs y
if plot_alpha
    for i_a = 1:length(i_aa)
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
            plt(i) = plot(wings_rect(i).LL_res.y/wings_rect(i).b*2, ...
                          wings_rect(i).LL_res.alpha_i(:,i_aa(i_a)), ...
                          LineWidth=lw(i), ...
                          Marker=markers(i), MarkerSize=ms(i), ...
                          DisplayName=disp_name);
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
            exp_name = fullfile(exp_fld, ...
                sprintf('T2_alpha_i_vs_y_AoA=%.2f.pdf', alpha(i_aa(i_a))));
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
    end
else
    disp('induced alpha vs y not plotted')
end

fig_count = fig_count + length(i_aa);

%Plot C_l vs alpha
if plot_C_l
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
        plt(i) = plot(wings_rect(i).LL_res.alpha, ...
                      wings_rect(i).LL_res.C_l_tot, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=disp_name);
    end
    plt(end+1) = plot(alpha, deg2rad(alpha-NACA_4415.alpha_L0)*NACA_4415.m_0, ...
                  LineWidth=1, color="k", LineStyle = "-.", ...
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
    ylabel('$C_l$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, 'T2_C_l_vs_alpha.pdf');
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
else
    disp('C_l vs alpha not plotted')
end

fig_count = fig_count + 1;

%Plot C_d vs alpha
if plot_C_di
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
        plt(i) = plot(wings_rect(i).LL_res.alpha, ...
                      wings_rect(i).LL_res.C_di_tot, ...
                      LineWidth=lw(i), ...
                      Marker=markers(i), MarkerSize=ms(i), ...
                      DisplayName=disp_name);
    end
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
    ylabel('$C_{d_i}$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
        if savefigs
            exp_name = fullfile(exp_fld, 'T2_C_di_vs_alpha.pdf');
            exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
                'BackgroundColor', 'none', 'Resolution', 300);
        end
else
    disp('C_{d_i}nduced vs alpha not plotted')
end

fig_count = fig_count + 1;