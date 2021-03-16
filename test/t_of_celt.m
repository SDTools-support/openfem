 % tests on of_celt (compiled m-type elements) :
 % bar1, beam1

 %------------------------------------------------- BAR1
 FEnode = [1  0 0 0  0 0 0;
           2  0 0 0  1 0 0 ];
 FEel0=[Inf abs('bar1');1 2 1001 1002 0];
 femesh('divide 4');
 model.Node=FEnode; model.Elt=FEel0;
 model.pl=m_elastic('dbval 1001 steel');
 model.il=p_beam('dbval 1002 rectangle .1 .1');

 [Case,model.DOF]=fe_mknl('init',model);
 k=fe_mknl('assemble',model,Case,1);

 %------------------------------------------------- BEAM1
 FEnode = [1  0 0 0  0 0 0;
           2  0 0 0  1 0 0 ];
 FEel0=[Inf abs('beam1');1 2 1001 1002 0 0];
 femesh('divide 4');
 model.Node=FEnode; model.Elt=FEel0;
 model.pl=m_elastic('dbval 1001 steel');
 model.il=p_beam('dbval 1002 rectangle .1 .1');

 [Case,model.DOF]=fe_mknl('init',model);
 k=fe_mknl('assemble',model,Case,1);

 %------------------------------------------------- BEAM1

 
