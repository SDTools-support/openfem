
if ~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt');return;end
wd=fullfile(fileparts(which('fe_mk')),'test');addpath(wd);

pw0=pwd;pw_=pw0;
diary off
delete(fullfile(tempdir,'openfemcheck.log'));
diary(fullfile(tempdir,'openfemcheck.log'));

t_mat_og
t_fe_curve

basic_elt_test('all')

basic_elt_test('compare')

% xxx
t_rotating
t_of_celt
t_fe_time 
t_thermal
t_quad4 
t_tria3
SL
