classdef EllipticWing
    % EllipticWing Represents an elliptic wing with specified aspect ratio and span.
    %   This class provides methods to calculate geometric properties of an elliptic wing,
    %   such as wing area, root chord, chord length, and spanwise coordinates.
    
    properties
        AR  % Aspect Ratio
        b   % Wing Span
        S   % Wing area
        c_root   % Root chord
        c_mean   % Mean chord
    end
    
    methods
        function obj = EllipticWing(AR, b)
            % EllipticWing Constructor
            %   obj = EllipticWing(AR, b) creates an elliptic wing object with the
            %   specified aspect ratio (AR) and wing span (b).
            %
            %   Inputs:
            %       AR - Aspect ratio (span^2 / wing area)
            %       b  - Wing span
            obj.AR = AR;
            obj.b = b;
            obj.S = obj.b^2 / obj.AR;
            obj.c_root = 4 * obj.S / (pi * obj.b);
            obj.c_mean = b/AR;
        end
        
        function [y, theta] = generate_coordinates(obj, N)
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
            theta = linspace(0, pi, N+2);
            theta = theta(2:end-1);  % Remove the first and last point to avoid singularities
            y = obj.b/2 * cos(theta);   % Spanwise locations
        end
        
        function c = chord_length(obj, y)
            % chord_length Calculate the chord length at a spanwise location.
            %   c = chord_length(obj, y) returns the chord length at the specified
            %   spanwise location (y) for an elliptic wing.
            %
            %   Input:
            %       y - Spanwise location
            %
            %   Output:
            %       c - Chord length at location y
            c = obj.c_root * sqrt(1 - (2 * y / obj.b).^2);
        end
    end
end

