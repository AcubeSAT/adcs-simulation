function [satrec, x]=orbit_init()
%% Get tle
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
infilename = "SSO-500-6PM.TLE";
%infilename = "SSO550_2.TLE";
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