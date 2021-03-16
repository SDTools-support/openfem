if ~(comstr(mexext,'dll')); return; end

rmpath('d:\dis_sdt\sdt.cur\sdtdev')
rmpath('d:\dis_sdt\sdt.cur\6.0')
rmpath('d:\dis_sdt\sdt.cur\openfem')
rmpath('d:\dis_sdt\sdt.cur')
rmpath('d:\dis_sdt\sdt.cur\sdtdemos')

a=getpref('OpenFEM','MexTarget');
%setpref('OpenFEM','MexTarget','d:\dis_sdt\sdt.cur\6.0')

test_dir='i:\tmp\openfem2\openfem';
setpref('OpenFEM','MexTarget',test_dir)

cd(test_dir)
ofutil('path')

