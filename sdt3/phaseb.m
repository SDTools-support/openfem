function phase=phaseb(xf,j1)

%PHASEB	Computes the phase of a complex matrix (in degrees)
%
%       Synopsis: phase=phaseb(xf)
%
%	The phase is unwrapped along each column unless a second arbitrary
%	argument is specified PHASE=PHASEB(XF,1).

%	Etienne Balmes  05/28/92, 03/30/94
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

[nw,ns] = size(xf);

phase=atan2(imag(xf),real(xf))/pi*180;
if nargin ==2; return; end

for j1 = 1:ns
  jumps = phase(2:nw,j1) - phase(1:nw-1,j1);
  ind   = [[find(abs(jumps)>330)]+1;nw+1];
  cor=0;

  for j2 = 1:(length(ind)-1)
	cor=cor-round(jumps(ind(j2)-1)/360)*360;
	if cor ~= 0
	 phase(ind(j2):ind(j2+1)-1,j1)=phase(ind(j2):ind(j2+1)-1,j1)+cor;
	end
  end
end

if nw>1; aa=mean(phase,1); else aa = phase; end
aa = aa-aa(1);
  for j1 = 2:ns
	if round(aa(j1)/360)~=0
	   phase(:,j1) = phase(:,j1) - round(aa(j1)/360)*360;
	end
  end
