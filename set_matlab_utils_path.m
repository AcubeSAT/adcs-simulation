function set_matlab_utils_path()

path = strrep(mfilename('fullpath'), 'set_matlab_utils_path','');

addpath([path 'src/utils/']);
addpath([path 'src/utils/quaternions/']);
addpath([path 'src/utils/system_model/']);
addpath([path 'src/utils/kf_lib/']);

import_kf_lib();

end
