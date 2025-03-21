classdef LiftingLine
    %LiftingLine Class to calulate the lift and drag of a wing using the
    %lifting line theory by Prandtl
    methods (Static)
        function A = solve_coeffs (wing, y, theta, alpha, N, m_0, alpha_L0)
            % solve_coeffs Solve for Fourier coefficients using Prandtl's Lifting Line Theory.
            %   A = solve_fourier_coefficients(wing, alpha, N) solves for the Fourier coefficients
            %   A_n using Prandtl's Lifting Line Theory for the specified wing, angle of attack,
            %   and number of spanwise points.
            %
            %   Inputs:
            %       wing - The wing to use for the calculation
            %       alpha - Angle of attack in degrees
            %       N - Number of terms in the Fourier series
            %       m_0 - Linear lift slope
            %       alpha_L0 - Zero lift angle of attack
            %
            %   Output:
            %   	y - Cartesian section position
            %       theta - Polar section positions
            %       A - Fourier coefficients
            if nargin<7
                alpha_L0 = -4;
            end
            alpha_L0 = deg2rad(alpha_L0);
        
            if nargin<6
                m_0 = 2*pi;
            end
            if nargin<5
                N = 50;
            end
        
            % Convert angle of attack to radians
            alpha = deg2rad(alpha);
            N_alpha = length(alpha);
        
            % Calculate chord lengths
            c_theta = wing.chord_length(y);
        
            % Initialize the matrix and right-hand side vector
            M = zeros(N, N);
            RHS = zeros(N, N_alpha);
            A = zeros(N, N_alpha);
        
            % Loop over each collocation point
            for i = 1:N
                % Loop over each Fourier coefficient
                for n = 1:N
                    M(i, n) = (4 * wing.b / (m_0 * c_theta(i))) * sin(n * theta(i)) ...
                               + n * sin(n * theta(i)) / sin(theta(i));
                end
                RHS(i, :) = alpha - alpha_L0;  % Angle of attack
            end
        
            % Solve for the Fourier coefficients
            for i = 1:N_alpha
                A(:, i) = M \ RHS(:, i);
            end
        end
        
        function [C_l, C_d] = calc_lift_drag_wing(wing, A)
            % calc_lift_drag Calculate the lift and drag coefficients of 
            % the entire wing based on the fourier coefficients.
            %
            %   Inputs:
            %       A - Fourier coefficients
            %
            %   Output:
            %       C_l - Lift coefficient
            %       C_d - Drag coefficient
        
            C_l = pi*wing.AR*A(1,:);
            C_d = pi*wing.AR*sum([1:size(A, 1)]'.*(A.^2), 1);
        end
        
        function [alpha_i, C_li, C_di, Gamma] = calc_lift_drag_sections(wing, y, theta, A)
            % calc_lift_drag Calculate the section lift and drag 
            % coefficients of based on the fourier coefficients.
            %
            %   Inputs:
            %       A - Fourier coefficients
            %
            %   Output:
            %       alpha - Angle of attack [deg]
            %       C_l - Lift coefficient
            %       C_d - Drag coefficient

            % Get airfoil chord lengths
            c_theta = wing.chord_length(y);

            % Vectorized calculation 
            % N = size(A, 1);
            % n_alpha = size(A,2);
            % n = repmat(reshape(1:N, [1, 1, N]), [N, n_alpha, 1]);
            % A_ext = repmat(reshape(A', 1, n_alpha, N), N, 1, 1);
            % theta_ext = repmat(theta', 1, n_alpha, N);
            % alpha_i = sum(n.*A_ext.*sin(n.*theta_ext)./sin(theta_ext), 3);
            % C_li = 4*wing.b./c_theta' .* sum(A_ext.*sin(n.*theta_ext), 3);
            % C_di = C_li .* alpha_i;
            % alpha_i = rad2deg(alpha_i);

            % Calculation with loops
            alpha_i = zeros(size(A));
            C_li = zeros(size(A));
            Gamma = zeros(size(A));
            N = size(A, 1);
            n_alpha = size(A,2);
            n = 1:N;
            for i_a = 1:n_alpha
                for i_t = 1:N
                    alpha_i(i_t, i_a) = sum(n.*A(:,i_a)'.*sin(n.*theta(i_t))./sin(theta(i_t)));
                    C_li(i_t, i_a) = 4*wing.b./c_theta(i_t) .* sum(A(:,i_a)'.*sin(n.*theta(i_t)));
                    Gamma(i_t, i_a) = 2*wing.b * sum(A(:,i_a)'.*sin(n.*theta(i_t)));
                end
            end
            C_di = C_li .* alpha_i;
            alpha_i = rad2deg(alpha_i);
        end
    end
end