%=========================================================%
%                      DEMO_NOPO                          %
%=========================================================%

%---------------------------------------------------------%
% 1. Importing models from other codes                    %
% See section 3.1.3 of the tutorial                       %
%---------------------------------------------------------%
pw0=pwd;
cd(fullfile(fileparts(which('fe_mk')),'demos'))
model = nopo('read -p 3d ex3d');
cd(pw0)
mflag=~system('medit');
if mflag medit('write ex3d',model);
else feplot(model);
end
