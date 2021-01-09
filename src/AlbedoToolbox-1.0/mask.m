% MASK Make masked array from two equal arrays.
%
% res = mask(refl_data,mask,contrast)
%
% refl_data is the reflectivity "background" data. mask is an array of
% logical values of visible indices. contrast [0..1] specifies the contrast of
% visible and not visible data points in refl_data.
%
% $Id: mask.m,v 1.4 2006/05/17 14:39:17 danji Exp $

function result = mask(refl_data,mask,contrast);

if nargin < 3
	contrast = 0.3;
end

result = (~mask*contrast+mask).*refl_data;

return
