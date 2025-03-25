classdef RectangularWing
    % RectangularWing Represents a rectangular wing with specified aspect ratio and span.
    %   This class provides methods to calculate geometric properties of a rectangular wing,
    %   such as wing area, chord length, and spanwise coordinates.
    
    properties
        AR  % Aspect Ratio
        b   % Wing Span
        S   % Wing area
        c_root   % Root chord
    end
    
    methods
        function obj = RectangularWing(AR, b)
            % RectangularWing Constructor
            %   obj = RectangularWing(AR, b) creates a rectangular wing object with the
            %   specified aspect ratio (AR) and wing span (b).
            %
            %   Inputs:
            %       AR - Aspect ratio (span^2 / wing area)
            %       b  - Wing span
            obj.AR = AR;
            obj.b = b;
            obj.S = obj.b^2 / obj.AR;
            obj.c_root = obj.S / obj.b;
        end
        
        function [y, theta] = generate_coordinates(obj, N, density_factor)
            % generate_coordinates Generate spanwise coordinates.
            %   [y, theta] = generate_coordinates(obj, N) generates the spanwise
            %   coordinates (y) and corresponding angular coordinates (theta) for
            %   N spanwise points.
            %
            %   Input:
            %       N - Number of spanwise points
            %
            %   Outputs:
            %       y     - Spanwise locations
            %       theta - Angular coordinates


            % theta = linspace(0, pi, N+2);
            % theta = theta(2:end-1);  % Remove the first and last point to avoid singularities
            %Calculate spanwise coordinates
            % y = -obj.b/2 * cos(theta);   % Spanwise locations

            if nargin<3
                density_factor = 0;
            end

            % density_factor = 1 + 1/(density_factor-.6)^2;
            if density_factor == 0
                theta = linspace(0, pi, N+2);
                theta = theta(2:end-1);  % Remove the first and last point to avoid singularities

                % Calculate spanwise coordinates
                y = -obj.b/2 * cos(theta);
            else
                x = linspace(-5, 5, N+2);
                y0 = (sqrt(4*density_factor.^2*x.^2 +1) - 1)./(2*x);
                y1 = y0./(max(abs(y0)));
                y1 = y1(2:end-1);
                y = y1*obj.b/2;
                theta = acos(-y1);
            end
        end
        
        function c = chord_length(obj, y)
            % chord_length Calculate the chord length at a spanwise location.
            %   c = chord_length(obj, y) returns the chord length at the specified
            %   spanwise location (y) for a rectangular wing.
            %
            %   Input:
            %       y - Spanwise location
            %
            %   Output:
            %       c - Chord length at location y
            c = obj.S / obj.b * ones(size(y));
        end
    end
end

