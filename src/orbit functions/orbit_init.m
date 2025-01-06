% -----------------------------------------------------------------------------
%
%                              procedure orbit_init
%
%  this procedure extracts all required data from the TLE file
%  and initializes the SGP4 orbit propagator
%
%   inputs        :
%   Asks for a TLE file
%
%   outputs       :
%     satrec      - struct including all required SPG4 orbit propagator variables
%     x           - random variable used if TLE random pick is enabled
%  ----------------------------------------------------------------------------*/


function [satrec, x] = orbit_init()

    %% Get tle

    %% Random TLE selection
    x = rand;
    % if x >= 0 && x <= 0.25
    %     infilename = "SSO-500-6PM.TLE";
    % elseif x > 0.25 && x <= 0.5
    %     infilename = "SSO-500-8PM.TLE";
    % elseif x > 0.5 && x <= 0.75
    %     infilename = "SSO-500-9PM.TLE";
    % else
    %     infilename = "SSO-500-11PM.TLE";
    % end

    %% Set TLE manually
    %infilename = "SSO-500-6PM.TLE";
    infilename = "SSO-600-11PM.TLE";

    %% Aks the user to provide TLE file
    %infilename = input('input elset filename: ','s');

    %% Open TLE
    infile = fopen(infilename, 'r');
    if (infile == -1)
        fprintf(1, 'Failed to open file: %s\n', infilename);
        return;
    end
    longstr1 = fgets(infile, 130);
    while ((longstr1(1) == '#') && (feof(infile) == 0))
        longstr1 = fgets(infile, 130);
    end
    longstr2 = fgets(infile, 130);

    %% Store tle values in satrec class and initialize orbit propagator
    satrec = twoline2rvMOD(longstr1, longstr2);

end
