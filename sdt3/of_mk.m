
% OF_MK OpenFEM gateway function for element matrix computations
%
% This function is only of interest for developpers
%
% Building of material caracteristic vector is done in fe_mat('of_mk')
%
% IOPT Integer options passed to of_mk follow the format
% [MatId EltId ndof Unused    100    Unused telt tcar NoRef] for loads
% [MatId EltId ndof SiteOfOut 200   Unused telt tcar]    for stress computation
% [MatId EltId ndof SiteOfOut MatOpt Unused telt tcar]   for matrix assembly
%
% telt : 1 plane stress, 2 plane strain, 3 axi
% tcar : 1 iso plane stress, 2 iso plane strain, 3 aniso
%
% CONSTIT follows the format
%  [Rho Eta StiffValues]
%
% For actual calling formats look at the source code of
%  fe_load, fe_stres, and element functions
%
%
% Compilation informations :
% 
% Choose compiler with : there are 2 choices
% mex -setup 
%
% Editing mexopt.sh is not necessary in most cases
%
% To compile with exact compiler options, prelink commands 
% (if appropriate), and linker options :
% mex -v toto.f
%       $Revision: 1.1 $  $Date: 2004/08/26 17:45:56 $
