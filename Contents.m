
% OpenFEM : an open Finite Element Library
% Version 2.0 beta2  Nov-2004
% Copyright (c) 2001-2002 Inria & SDTools
%
% OpenFEM is a registered trademark. You can use it under the terms given
% in trademark.html in this folder or http://www-rocq.inria.fr/OpenFEM
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License version 2.1 as published by the Free Software Foundation
% See LGPL.txt in this folder or http://www.gnu.org/copyleft/lesser.html
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% For more details about the library or to contact the authors
%
%   www.inria.fr/openfem
%   www.sdtools.com/openfem
%
%
% 2-D Elements
%
%   q4p        - 4 node 8 DOF quadrangle
%   q5p        - 5 node 10 DOF quadrangle      
%   q8p        - 8 node 16 DOF quadrangle      
%   t3p        - 3 node 6 DOF triangle      
%   t6p        - 6 node 12 DOF triangle
%
% 3-D plate/shell
%
%   dktp       - 3 node 9 DOF discrete Kirchoff plate
%   mitc4      - 4 node 20 DOF shell
%   quad4      - quadrilateral 4-node 20/24-DOF plate/shell element
%   quad9      - (display only)
%   quadb      - quadrilateral 8-node 40/48-DOF plate/shell element
%   tria3      - 3-node 15/18-DOF thin plate/shell element
%   tria6      - 6-node 30/36-DOF plate/shell element (display only)
%  
% 3-D Volumes
%
%   hexa27     - 27-node 81-DOF isoparametric solid element
%   hexa20     - 20-node 60-DOF isoparametric solid element
%   hexa8      - 8-node 24-DOF isoparametric solid element
%   penta15    - 15 node 45-DOF solid element
%   penta6     - 6-node 18-DOF solid element
%   tetra10    - 10-node 30-DOF solid element
%   tetra4     - 4-node 12-DOF isoparametric solid element
%
% 3-D Volumes (supporting non linearities)
%
%   hexa27b    - 27-node 81-DOF isoparametric solid element
%   hexa20b    - 20-node 60-DOF isoparametric solid element
%   hexa8b     - 8-node 24-DOF isoparametric solid element
%   penta15b   - 15 node 45-DOF solid element
%   penta6b    - 6-node 18-DOF solid element
%   tetra10b   - 10-node 30-DOF solid element
%   tetra4b    - 4-node 12-DOF isoparametric solid element
%
% Other elements
%
%   bar1       - standard 2-node 6-DOF bar element
%   beam1      - standard 2-node 12-DOF Bernoulli-Euler beam element
%   beam1      - pretensionned 2-node 12-DOF Bernoulli-Euler beam element
%   beam3      - display only
%   celas      - scalar springs and penalized rigid links
%   mass1      - concentrated mass/inertia element
%   mass2      - concentrated mass/inertia element with offsets
%   rigid      - handling of linearized rigid links
%
% Miscellaneous utilities
%
%   basis      - coordinate transformation utilities
%   getegroup  - find group headers in a model description matrix
%   ofutil     - developper utilities for the OpenFEM library
%   of_mk      - Interface DLL for FORTRAN routines
%   readnopo   - reading of MODULEF NOPO format
%   stack_get, stack_rm, stack_set - stack handling utilities
%   sdtw       - extended warning handling
%
% General FEM utilities
%
%   fe_c          - DOF selection and I/O shape matrix construction
%   fe_load       - creation of distributed load vectors
%   fe_mat        - material property handling
%   fe_mk         - assembly of full and reduced FE model matrices
%   fe_mkn;       - assembly optimized for non-linear operation
%   fe_stres      - element energy and stress computations
%   feplot,fecom  - GUI commands for manipulation of deformation plots
%   femesh,feutil - UI commands for finite element mesh handling
%   feplot        - GUI for deformed structure visualization
%   feutil        - mesh handling utilities
