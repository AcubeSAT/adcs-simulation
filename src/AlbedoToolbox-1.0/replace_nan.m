% REPLACE_NAN Replaces NaN values with annual mean or specified
% reflectivity data.
%
% new_refl = replace_nan(main_refl [, param])
%
% main_refl is the reflectivity data in which undefined values will be
% replaced. If param is a a REFL struct, data from param is copied to
% main_refl. If param is a string it is used as a path to the library of
% annual reflectivity data. The libaray should contain a directory for each
% year of data, and the files should must be of the structure'gaYYMMDD.mat'
% and the annual mean data should be 'gaYY0101-YY1231.mat'. Per default,
% i.e. only main_refl is specified, the path is searched for the annual
% mean.
%
% $Id: replace_nan.m,v 1.6 2006/05/17 14:39:17 danji Exp $

function new_refl = replace_nan(main_refl, param);

persistent nan_refl;

% Check parameters
if nargin > 1 && isstruct(param)
    index = find(isnan(main_refl.data));
    main_refl.data(index) = param.data(index);
    msg = 'parameter.';
else
  year = jd2d(main_refl.start_time);
  if ~(isfield(nan_refl,'data') && year == jd2d(nan_refl.start_time))
    if year > 2000
      fyear = year - 2000;
    else
      fyear = year - 1900;
    end
    filename = sprintf('ga%02.0f0101-%02.0f1231.mat',fyear,fyear);
    file = which(filename);
    if isempty(file) && exist('param')
      filename = fullfile(param,num2str(year),filename);
      msg = filename;
    else
      msg = file;
    end
    nan_refl = load(filename);
  else
    msg = 'pre-loaded data.';
  end
  index = find(isnan(main_refl.data));
  main_refl.data(index) = nan_refl.data(index);
end

new_refl = main_refl;
if size(index,1) > 0
  disp(['replace_nan.m: ' num2str(size(index,1)) ' NaN replacements from ' msg]);
  new_refl.type = 'Mean supported raw';
end

return
