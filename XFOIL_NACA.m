% This class provides the option to run XFOIL from Matlab for a 
% specific airfoil and read the results form the output files.
classdef XFOIL_NACA
    properties
        NACA
        xfoil_exe
        Re
        Ma
        n_crit
        Xt
    end
    methods
        %% --- Class constructor
        function obj = XFOIL_NACA(NACA, xfoil_exe, Re, Ma, n_crit, Xt)
            obj.NACA = NACA;
            if nargin<2
                xfoil_exe = fullfile('..','..', ...
                    '03_additional_material', 'Xfoil', 'xfoil.exe');  %Path to Xfoil executable
            end
            obj.xfoil_exe = xfoil_exe;

            if nargin<3
                Re = 6e6;
            end
            obj.Re = Re;

            if nargin<4
                Ma = 0;
            end
            obj.Ma = Ma;

            if nargin<5
                n_crit = 9;
            end
            obj.n_crit = n_crit;

            if nargin<6
                Xt = .1;
            end
            obj.Xt = Xt;

        end
        
        %% --- Function to calculate C_l and C_d of 4-digit NACA for a range of angles of attack---
        function [alpha, C_l, C_d, xt_top, xt_bot] = calc_c_ld (obj, free_BL, use_cached_input, ...
                AoAs, numNodes, max_iter)
            if nargin < 6
                max_iter = 200;
            end
            if nargin < 5
                numNodes = 160;
            end
            if nargin < 4
                AoAs = [-5, 15, .25];  % Angles of attack [deg]
            end
            if nargin < 3
                use_cached_input = false;
            end
            if nargin < 2
                free_BL = false;
            end
            
            % Xfoil preparation & calculation
            exp_fld = 'xfoil_exports' 
            if free_BL
                fname_input = ['input_NACA' obj.NACA '_Cl_free_BL.txt'];  % XFoil input filename
                fname_res = ['Cl_NACA' obj.NACA '_free_BL.txt'];  % Lift coefficient filename
            else
                fname_input = ['input_NACA' obj.NACA '_Cl_fixed_BL.txt'];  % Airfoil coordinates filename
                fname_res = ['Cl_NACA' obj.NACA '_fixed_BL.txt'];  % Lift coefficient filename
            end
            
            fname_coords = ['Coords_NACA' obj.NACA '.txt'];  % Airfoil coordinates filename
            
            % Create export directory if it doesn't exist
            if ~exist(exp_fld, 'dir')
                mkdir(exp_fld);
            end

            % Delete files if they exist
            if (exist(fullfile(exp_fld, fname_coords),'file'))
                delete(fullfile(exp_fld, fname_coords));
            end
            if (exist(fullfile(exp_fld, fname_res),'file'))
                delete(fullfile(exp_fld, fname_res));
            end

            if use_cached_input
                assert(exist(fullfile(exp_fld, fname_input),'file'), "Input file not found")
            else
                % Prepare input to xfoil
                % Create the airfoil
                fid = fopen(fullfile(exp_fld, fname_input), 'w');
                fprintf(fid, ['NACA ' obj.NACA '\t\t\t\t\t! Airfoil selection\n']);
                fprintf(fid, 'PPAR\t\t\t\t\t\t! show paneling menu\n');
                fprintf(fid, sprintf('N %d\t\t\t\t\t\t! Number of panels\n', numNodes));
                fprintf(fid, '\n\n');
                
                % Save the airfoil data points
                fprintf(fid,['PSAV ' exp_fld '\\' fname_coords '\n']);
                
                % Calculate C_l vs alpha values
                fprintf(fid, 'OPER\t\t\t\t\t\t! Enter operational menu\n');
                fprintf(fid, sprintf('ITER %d\t\t\t\t\t\t\t! Maximum no. of iterations\n', max_iter));
                fprintf(fid, sprintf('VISC %10e\t\t\t\t! Set Reynolds number\n', obj.Re));
                fprintf(fid, sprintf('Mach %.4f\t\t\t\t\t! Set Mach number\n', obj.Ma));
                
                fprintf(fid,'VPAR\t\t\t\t\t\t! Enter BL parameter menu\n');
                if free_BL
                    fprintf(fid, sprintf(['N %d\t\t\t\t\t\t\t! '...
                        'Set critical amplification exponent\n'], obj.n_crit));
                else
                    fprintf(fid, sprintf(['N %d\t\t\t\t\t\t\t! '...
                        'Set critical amplification exponent\n'], obj.n_crit));
                    fprintf(fid, 'XTR\n');
                    fprintf(fid, sprintf('%.3f\n', obj.Xt));  % Set upper trip positions
                    fprintf(fid, '\n');
                end
                fprintf(fid, '\n');
                
                % Configure output file
                fprintf(fid, 'PACC\t\t\t\t\t\t! Enable polar accumulation (stores results)\n');
                fprintf(fid, [exp_fld '\\' fname_res '\n']);
                fprintf(fid, '\n');
                
                %Run calculation for all AoAs
                fprintf(fid, sprintf('ASEQ %s\n' , sprintf('%.2f ' , AoAs)));
                
                %Finish up
                fprintf(fid, 'PACC\t\t\t\t\t\t! Disable polar accumulation (closes output file)\n');
                fprintf(fid, '\nQUIT\n');
                
                % Close file
                fclose(fid);

                % Run XFoil using input file
                cmd = obj.xfoil_exe + " < " + fullfile(exp_fld, fname_input);
                %[status,result] = system(cmd);
                system(cmd)

                [alpha, C_l, C_d, xt_top, xt_bot] = obj.read_C_ld (free_BL);
            end
        end

        %% --- Function to calculate C_p of 4-digit NACA---
        function [x, y, C_p, x_trans] = calc_c_p (obj, free_BL, use_cached_input, ...
                AoA, numNodes, max_iter)
            if nargin < 6
                max_iter = 200;
            end
            if nargin < 5
                numNodes = 500;
            end
            if nargin < 4
                AoA = 10;  % Angle of attack [deg]
            end
            if nargin < 3
                use_cached_input = false;
            end
            if nargin < 2
                free_BL = false;
            end
            
            % Xfoil preparation & calculation
            exp_fld = 'xfoil_exports' 
            if free_BL
                fname_input = ['input_NACA' obj.NACA '_Cp_free_BL.txt'];  % XFoil input filename
                fname_res = ['Cp_NACA' obj.NACA '_free_BL.txt'];  % Lift coefficient filename
                fname_trans = ['trans_NACA' obj.NACA '_free_BL.txt'];  % Transition point filename
            else
                fname_input = ['input_NACA' obj.NACA '_Cp_fixed_BL.txt'];  % Airfoil coordinates filename
                fname_res = ['Cp_NACA' obj.NACA '_fixed_BL.txt'];  % Lift coefficient filename
            end
            
            fname_coords = ['Coords_NACA' obj.NACA '.txt'];  % Airfoil coordinates filename
            
            % Create export directory if it doesn't exist
            if ~exist(exp_fld, 'dir')
                mkdir(exp_fld);
            end

            % Delete files if they exist
            if (exist(fullfile(exp_fld, fname_coords),'file'))
                delete(fullfile(exp_fld, fname_coords));
            end
            if (exist(fullfile(exp_fld, fname_res),'file'))
                delete(fullfile(exp_fld, fname_res));
            end
            if free_BL
                if (exist(fullfile(exp_fld, fname_trans),'file'))
                delete(fullfile(exp_fld, fname_trans));
            end
            end

            if use_cached_input
                assert(exist(fullfile(exp_fld, fname_input),'file'), "Input file not found")
            else
                % Prepare input to xfoil
                % Create the airfoil
                fid = fopen(fullfile(exp_fld, fname_input), 'w');
                fprintf(fid, ['NACA ' obj.NACA '\t\t\t\t\t! Airfoil selection\n']);
                fprintf(fid, 'PPAR\t\t\t\t\t\t! show paneling menu\n');
                fprintf(fid, sprintf('N %d\t\t\t\t\t\t! Number of panels\n', numNodes));
                fprintf(fid, '\n\n');
                
                % Save the airfoil data points
                fprintf(fid,['PSAV ' exp_fld '\\' fname_coords '\n']);
                
                % Calculate C_l vs alpha values
                fprintf(fid, 'OPER\t\t\t\t\t\t! Enter operational menu\n');
                fprintf(fid, sprintf('ITER %d\t\t\t\t\t\t\t! Maximum no. of iterations\n', max_iter));
                fprintf(fid, sprintf('VISC %10e\t\t\t\t! Set Reynolds number\n', obj.Re));
                fprintf(fid, sprintf('Mach %.4f\t\t\t\t\t! Set Mach number\n', obj.Ma));
                
                fprintf(fid,'VPAR\t\t\t\t\t\t! Enter BL parameter menu\n');
                if free_BL
                    fprintf(fid, sprintf(['N %d\t\t\t\t\t\t\t! '...
                        'Set critical amplification exponent (free transition)\n'], obj.n_crit));
                else
                    fprintf(fid, 'XTR\n');
                    fprintf(fid, sprintf('%.3f\n', obj.Xt));  % Set upper trip positions
                    fprintf(fid, '\n');
                end
                fprintf(fid, '\n');
                
                %Run calculation for 
                fprintf(fid, sprintf('ALFA %.3f\n' , AoA));

                % Save output to txt file
                fprintf(fid, ['CPWR ' exp_fld '\\' fname_res '\n']);
                if free_BL
                    fprintf(fid, 'vplo\n');  % Enter BL plot menu
                    fprintf(fid, 'N\n');  % Enter BL plot menue
                    fprintf(fid, ['dump ' exp_fld '\\' fname_trans '\n']);  % Save transition data
                    fprintf(fid,'\n');
                end
                fprintf(fid,'\n');
                
                %Finish up
                fprintf(fid, '\nQUIT\n');
                
                % Close file
                fclose(fid);

                % Run XFoil using input file
                cmd = obj.xfoil_exe + " < " + fullfile(exp_fld, fname_input);
                %[status,result] = system(cmd);
                system(cmd)

                [x, y, C_p, x_trans] = obj.read_C_p (free_BL);
            end
        end

        %% --- Function to calculate ΔC_p of 4-digit NACA---
        function [xc, dC_p] = calc_dc_p (obj, x, y, C_p)
            % Split upper and lower surfaces
            upper_idx = y >= 0;
            lower_idx = y < 0;

            % Sort upper and lower surfaces by x-coordinates
            [x_upper, sort_idx_upper] = sort(x(upper_idx));
            C_p_upper = C_p(upper_idx);
            C_p_upper = C_p_upper(sort_idx_upper);
            
            [x_lower, sort_idx_lower] = sort(x(lower_idx));
            C_p_lower = C_p(lower_idx);
            C_p_lower = C_p_lower(sort_idx_lower);
            
            % Ensure x-coordinates match for interpolation
            C_p_lower_interp = interp1(x_lower, C_p_lower, x_upper, ...
                'linear', 'extrap');
            
            % Compute relative pressure coefficient (ΔCp)
            dC_p = C_p_lower_interp - C_p_upper;
            xc = x_upper;
        end

        %% --- Function to load the xy coordinates of a 4-digit NACA from a text file from Xfoil---
        function [XB, YB] = read_coords (obj)
            % Determine filename
            fpath_coords = ['xfoil_exports\Coords_NACA' obj.NACA '.txt'];

            % Read airfoil coordinates
            coords_airfoil_file = fopen(fpath_coords);

            dataBuffer = textscan(coords_airfoil_file,'%f %f','CollectOutput',1,...
                                             'Delimiter','','HeaderLines',0);
            fclose(coords_airfoil_file);

            % Separate boundary points
            XB = dataBuffer{1}(:,1);
            YB = dataBuffer{1}(:,2); 
        end

        %% --- Function to load the AoAs, C_l and C_d values of a 4-digit NACA from a text file from Xfoil---
        function [alpha, C_l, C_d, xt_top, xt_bot] = read_C_ld (obj, free_BL)
            % Determine filename
            if free_BL
                fpath_res = ['xfoil_exports\Cl_NACA' obj.NACA '_free_BL.txt'];
            else
                fpath_res = ['xfoil_exports\Cl_NACA' obj.NACA '_fixed_BL.txt'];
            end

            % Read lift coefficients
            res_file = fopen(fpath_res);
            dataBuffer = textscan(res_file,'%f %f %f %f %f %f %f', ...
                                        'HeaderLines',12,...
                                        'CollectOutput',1,...
                                        'Delimiter','');
            fclose(res_file);
            
            % Separate Cp data
            alpha  = dataBuffer{1,1}(:,1); 
            C_l  = dataBuffer{1,1}(:,2); 
            C_d = dataBuffer{1,1}(:,3);
            xt_top = dataBuffer{1,1}(:,6);
            xt_bot = dataBuffer{1,1}(:,7);
        end

        %% --- Function to load the x, y, and C_p values of a 4-digit NACA from a text file from Xfoil---
        function [x, y, C_p, x_trans] = read_C_p (obj, free_BL)
            % Determine filename
            if free_BL
                fpath_res = ['xfoil_exports\Cp_NACA' obj.NACA '_free_BL.txt'];
            else
                fpath_res = ['xfoil_exports\Cp_NACA' obj.NACA '_fixed_BL.txt'];
            end

            % Read lift coefficients
            res_file = fopen(fpath_res);
            dataBuffer = textscan(res_file,'%f %f %f', ...
                                        'HeaderLines',3,...
                                        'CollectOutput',1,...
                                        'Delimiter','');
            fclose(res_file);

            % Separate Cp data
            x  = dataBuffer{1,1}(:,1); 
            y  = dataBuffer{1,1}(:,2); 
            C_p = dataBuffer{1,1}(:,3);

            % Find transition point for free transition
            if free_BL
                res_BL_file = fopen(['xfoil_exports\trans_NACA' obj.NACA '_free_BL.txt']);
                dataBuffer = textscan(res_BL_file,'%f %f', ...
                                        'HeaderLines',7,...
                                        'CollectOutput',1,...
                                        'Delimiter','');
                fclose(res_BL_file);

                %Find transition point
                i_zero = find(dataBuffer{1,1}(:,2) == 0,2, 'first'); % find the first two points where nc is zero
                i_trans = i_zero(2)-1; % Index of the transition point on the upper surface
                x_trans = dataBuffer{1,1}(i_trans,1);
            else
                x_trans = .1;
            end    
        end
        
        %% --- Function to plot the shape of an airfoil---
        function plot_airfoil (obj, XB, YB)
            % Split airfoil into (U)pper and (L)ower
            XB_U = XB(YB >= 0);
            XB_L = XB(YB < 0);
            YB_U = YB(YB >= 0);
            YB_L = YB(YB < 0);
            
            % Plot Airfoil
            figure(1);
            cla; hold on; grid off;
            set(gcf,'Color','White');
            set(gca,'FontSize',12);
            plot(XB_U,YB_U,'b.-');
            plot(XB_L,YB_L,'r.-');
            xlabel('X Coordinate');
            ylabel('Y Coordinate');
            axis equal;
        end
        
        %% --- Function to plot the lift coefficient of an airfoil over a range of AoAs---
        function plot_Cl (obj, alpha, C_l)
            figure(2);
            cla; hold on; grid on;
            set(gcf,'Color','White');
            set(gca,'FontSize',12);
            
            plot(alpha , C_l, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
            
            % Highlight x=0 and y=0 grid lines
            gray_color = [0.2, 0.2, 0.2];
            xline(0, 'Color', gray_color, 'LineWidth', 1.5); % Thick vertical line at x=0
            yline(0, 'Color', gray_color, 'LineWidth', 1.5); % Thick horizontal line at y=0
            
            % Plot labels
            xlabel('AoA [deg]', 'Interpreter', 'latex');
            ylabel('$C_l$', 'Interpreter', 'latex');
            set(gca, 'TickLabelInterpreter', 'latex');
            
            ylim('auto');
            xticks(-10:2:15);
        end

        %% --- Function to plot the pressure coefficient of an airfoil over x ---
        function plot_Cp (obj, x, C_p)
            figure(3);
            cla; hold on; grid on;
            set(gcf,'Color','White');
            set(gca,'FontSize',12);
            
            plot(alpha , C_p, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
            
            % Highlight x=0 and y=0 grid lines
            gray_color = [0.2, 0.2, 0.2];
            xline(0, 'Color', gray_color, 'LineWidth', 1.5); % Thick vertical line at x=0
            yline(0, 'Color', gray_color, 'LineWidth', 1.5); % Thick horizontal line at y=0
            
            % Plot labels
            xlabel('$x/c$', 'Interpreter', 'latex');
            ylabel('$C_p$', 'Interpreter', 'latex');
            set(gca, 'TickLabelInterpreter', 'latex');
            
            ylim('auto');
            xticks(-1:1:.1);
        end
    end
end