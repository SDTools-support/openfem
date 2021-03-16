function fegui(in1)

%FEGUI	Verification that global variables for FEMESH are properly declared
%
%	FEMESH uses a number of global variables
%
%	  FEnode/FEn0/FEn1    main/selected/alternate node set
%	  FEelt/FEel0/FEel1   main/selected/alternate model description matrix
%
%       FEGUI checks that these variables are properly declared as global in
%	your base and caller function workspaces. 
%
%       USER CALLS TO FEGUI ARE OBSOLETE. USE FEMESH INSTEAD.
%
%	See also help femesh, feplot, fecom

%       Etienne Balmes  02/02/94, 07/26/02
%       Copyright (c) 1990-2020 by SDTools
%       All Rights Reserved.
%       $Revision: 1.3 $  $Date: 2020/02/26 08:47:47 $


evalin('caller','iigui({''FEnode'',''FEn0'',''FEn1'',''FEelt'',''FEel0'',''FEel1''});');
