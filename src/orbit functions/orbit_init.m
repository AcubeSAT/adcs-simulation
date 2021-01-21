function [satrec]=orbit_init()
%% Get tle
infilename = "SSO-500-6PM.TLE";
%infilename = "SSO500_29_08.TLE";
%infilename = input('input elset filename: ','s');
infile = fopen(infilename, 'r');
if (infile == -1)
        fprintf(1,'Failed to open file: %s\n', infilename);
        return;
end
longstr1 = fgets(infile, 130);
while ( (longstr1(1) == '#') && (feof(infile) == 0) )
    longstr1 = fgets(infile, 130);
end       
longstr2 = fgets(infile, 130);

%% Store tle values in satrec class
satrec = twoline2rvMOD(longstr1,longstr2);
end