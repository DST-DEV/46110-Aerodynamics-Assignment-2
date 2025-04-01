% This script calculates the lift slope and zero lift angle using XFOIL.
clear; clc; close all;
xfoil_exe = fullfile('..','..', ...
                    '03_additional_material', 'Xfoil', 'xfoil.exe');  %Path to Xfoil executable

%% Run Calculations
rerun_XFOIL = false;
plot_XFOIL = true

AoAs = [-20, 20, .05];
numNodes = 160;
max_iter = 400;

% If output file already exists, load it
exp_fld = 'results';
fpath_res = fullfile(exp_fld, 'XFOIL_NACA_4415.mat');
if ~exist(exp_fld, 'dir')
    mkdir(exp_fld);
elseif exist(fpath_res,'file')
    NACA_4415 = load(fpath_res).NACA_4415;
else
    NACA_4415 = struct();
end


%% Calculate C_l
if rerun_XFOIL
    XFOIL_4415 = XFOIL_NACA('4415', xfoil_exe, 6e6, 0, 9, .1);
    [alpha, C_l, C_d, xt_top, xt_bot] = ...
        XFOIL_4415.calc_c_ld(true, false, AoAs, numNodes, max_iter);

    NACA_4415.xfoil_res = struct('alpha', alpha, 'C_l', C_l, ...
            'C_d', C_d, 'xt_top', xt_top, 'xt_bot', xt_bot);

    save(fpath_res, 'NACA_4415');  % Save results
else
    alpha = NACA_4415.xfoil_res.alpha;
    C_l = NACA_4415.xfoil_res.C_l;
    C_d = NACA_4415.xfoil_res.C_d;
end

%% Calculate zero lift angle and lift slope
% Fit a linear model to the data (CL = p(1) * alpha + p(2))
i_sel = find(alpha>-10 & alpha<5);
alpha_sel = alpha(i_sel);
p = polyfit(alpha_sel, C_l(i_sel), 1);

% Get lift slope in [1/rad] (first coefficient of the linear polynomial) 
NACA_4415.m_0 = p(1) / (pi / 180);
% Get zero lift angle in [deg](root of the polynomial => y-intercept / slope
NACA_4415.alpha_L0 = -p(2) / p(1);
save(fpath_res, 'NACA_4415');  % Save results

if plot_XFOIL
    % Plot C_l
    figure(1);
    plot(alpha, C_l, 'b', 'DisplayName', 'Data'); % Plot data points

    hold on;
    alpha_fit = linspace(min(alpha_sel), max(alpha_sel), 50);
    C_l_fit = polyval(p, alpha_fit);
    plot(alpha_fit, C_l_fit, 'r-', 'DisplayName', 'Fitted Line');
    hold off
    
    xlabel('AoA [deg]');
    ylabel('C_l');
    legend;
    grid on;

    % Plot C_d
    figure(2);
    cla; hold on; grid on;
    plot(alpha, C_d, 'b', 'DisplayName', 'Data'); % Plot data points
    hold off;

    xlabel('AoA [deg]');
    ylabel('C_d');
    legend;
    grid on;
end