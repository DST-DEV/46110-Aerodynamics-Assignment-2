clear; clc;
%% Load results
res_ell = load('results\T1_wings_ell.mat').wings_ell;
res_rect = load('results\T2_wings_rect.mat').wings_rect;
res_taper = load('results\T3_wings_tapered.mat').wings_taper;
res_twist = load('results\T4_wings_rect_twist.mat').wings_twist;

%% Plots
rho = 1.225;
U_0 = 1;

% Selection
savefigs = true;
plot_Gamma = true;
plot_C_l = true;


% Settings
cols = ["#0072BD", "#D95319", "#EDB120", "#77AC30", "#80B3FF"];  % Colors of the lines
markers = ["none", "none", "none", "none", "none"];  % Markers for the four methods
ms = [4.5, 4.5, 4.5, 4.5, 4.5];  % Marker size for the plots of the four methods
lw = [1.5, 1.5, 1.5, 1.5, 1.5];  % Linewidth for the lines of the four methods
ax_col = [0.2, 0.2, 0.2];  % Color of accented axes
ax_lw = 1.5;  % Line width of accented axes
fs = 16;  % Plot font size
fig_count = 0;

% Create export directory if it doesn't exist
exp_fld = 'plots';
if ~exist(exp_fld, 'dir')
    mkdir(exp_fld);
end

%Plot Gamma non-dimensionalized vs y for AR = 6
if plot_Gamma
    % Create plot
    figure(fig_count+1);
    set(gcf, 'Position', [100, 100, 800, 400]);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot Tapered wing
    for i = 1:numel(res_taper)
        plot(res_taper(i).LL_res.y/res_taper(i).wing.b*2, ...
             res_taper(i).LL_res.Gamma_nd, ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=sprintf('Tapered wing $TR=%.1f$', res_taper(i).TR));
    end

    % Plot elliptic wing
    i_AR_ell = find([res_ell.AR] == 6);
    i_aoa_ell = find(res_ell(i_AR_ell).LL_res.alpha == 6);
    plot(res_ell(i_AR_ell).LL_res.y/res_ell(i_AR_ell).wing.b*2, ...
         res_ell(i_AR_ell).LL_res.Gamma_nd(:, i_aoa_ell), ...
         LineWidth=1.5, LineStyle = '--', color = 'k', ...
         DisplayName="Elliptic wing");

    % Plot rectangular wing
    i_AR_rect = find([res_rect.AR] == 6);
    i_aoa_rect = find(res_rect(i_AR_rect).LL_res.alpha == 6);
    plot(res_rect(i_AR_rect).LL_res.y/res_rect(i_AR_rect).wing.b*2, ...
         res_rect(i_AR_rect).LL_res.Gamma_nd(:, i_aoa_rect), ...
         LineWidth=1.5, LineStyle = ':', color = 'k', ...
         DisplayName="Rectangular wing");
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'eastoutside', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$\Gamma/(U_\infty\bar{c})$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'Comparison_T1-3_Gamma_nd_vs_y.pdf');
        exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('Gamma vs y not plotted');
end

fig_count = fig_count +1;

%Plot C_l vs y for AR = 6
if plot_C_l
    % Create plot
    figure(fig_count+1);
    set(gcf, 'Position', [100, 100, 800, 400]);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot Tapered wing
    for i = 1:numel(res_taper)
        plot(res_taper(i).LL_res.y/res_taper(i).wing.b*2, ...
             res_taper(i).LL_res.C_l, ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=sprintf('Tapered wing $TR=%.1f$', res_taper(i).TR));
    end

    % Plot elliptic wing
    i_AR_ell = find([res_ell.AR] == 6);
    i_aoa_ell = find(res_ell(i_AR_ell).LL_res.alpha == 6);
    plot(res_ell(i_AR_ell).LL_res.y/res_ell(i_AR_ell).wing.b*2, ...
         res_ell(i_AR_ell).LL_res.C_l(:, i_aoa_ell), ...
         LineWidth=1.5, LineStyle = '--', color = 'k', ...
         DisplayName="Elliptic wing");

    % Plot rectangular wing
    i_AR_rect = find([res_rect.AR] == 6);
    i_aoa_rect = find(res_rect(i_AR_rect).LL_res.alpha == 6);
    plot(res_rect(i_AR_rect).LL_res.y/res_rect(i_AR_rect).wing.b*2, ...
         res_rect(i_AR_rect).LL_res.C_l(:, i_aoa_rect), ...
         LineWidth=1.5, LineStyle = ':', color = 'k', ...
         DisplayName="Rectangular wing");
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'eastoutside', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$C_l$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'Comparison_T1-3_C_l_vs_y.pdf');
        exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('C_l vs y not plotted');
end

fig_count = fig_count +1;


%Plot C_l*c vs y for AR = 6
if plot_C_l
    % Create plot
    figure(fig_count+1);
    set(gcf, 'Position', [100, 100, 800, 400]);
    cla; hold on; grid on;
    colororder(cols);
    ax = gca;
        
    % Highlight y=0 grid line
    y_ax = xline(0, Color=ax_col, LineWidth=ax_lw, ...
                 HandleVisibility='off'); % Thick vertical line at x=0

    % Plot Tapered wing
    for i = 1:numel(res_taper)
        plot(res_taper(i).LL_res.y/res_taper(i).wing.b*2, ...
             res_taper(i).LL_res.C_l.*res_taper(i).wing.chord_length(res_taper(i).LL_res.y)', ...
             LineWidth=lw(i), Marker=markers(i), MarkerSize=ms(i), ...
             DisplayName=sprintf('Tapered wing $TR=%.1f$', res_taper(i).TR));
    end

    % Plot elliptic wing
    i_AR_ell = find([res_ell.AR] == 6);
    i_aoa_ell = find(res_ell(i_AR_ell).LL_res.alpha == 6);
    plot(res_ell(i_AR_ell).LL_res.y/res_ell(i_AR_ell).wing.b*2, ...
         res_ell(i_AR_ell).LL_res.C_l(:, i_aoa_ell).*res_ell(i_AR_ell).wing.chord_length(res_ell(i_AR_ell).LL_res.y)', ...
         LineWidth=1.5, LineStyle = '--', color = 'k', ...
         DisplayName="Elliptic wing");

    % Plot rectangular wing
    i_AR_rect = find([res_rect.AR] == 6);
    i_aoa_rect = find(res_rect(i_AR_rect).LL_res.alpha == 6);
    plot(res_rect(i_AR_rect).LL_res.y/res_rect(i_AR_rect).wing.b*2, ...
         res_rect(i_AR_rect).LL_res.C_l(:, i_aoa_rect).*res_rect(i_AR_rect).wing.chord_length(res_rect(i_AR_rect).LL_res.y)', ...
         LineWidth=1.5, LineStyle = ':', color = 'k', ...
         DisplayName="Rectangular wing");
    hold off; 

    % Configure limits and ticks
    ylim('auto');
    xticks(-1:.5:1);
    xlim(ax, [-1, 1]);

    % Plot labels
    set(gcf,'Color','White');
    set(ax,'FontSize',fs);
    legend('Location', 'eastoutside', 'Interpreter', 'latex')
    xlabel('$y/(b/2)$', 'Interpreter', 'latex');
    ylabel('$C_l$', 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');

    % Save figure
    if savefigs
        exp_name = fullfile(exp_fld, 'Comparison_T1-3_C_l_vs_y.pdf');
        exportgraphics(gcf, exp_name, 'ContentType', 'vector', ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
else
    disp('C_l vs y not plotted');
end

fig_count = fig_count +1;