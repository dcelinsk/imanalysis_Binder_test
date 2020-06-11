function positionVec=decodeShaftEncoder(encVlt,inchRadius)

% pass a rotary encoder's (not quadrature) values
% and get back a positionVec like a quadradture would.
% i ask for radius of load (wheel) in inches to give 
% back position in meters.

rInMet=inchRadius*0.0254;
wC=2*pi*rInMet;


% get range
%dR=max(encVlt)-min(encVlt);
%prefer this, but not essential  FROM ERIC
dR=max(encVlt);

% min meter delta
metPerInc=wC/dR;

% Make indxing easier by adding the zero here
difTM=[0,diff(encVlt)];

% roll overs.
%pRolInds=find(diff(encVlt)>(dR/2));   FROM ERIC
%I just like this indexing better , but not essential  FROM ERIC
fJumps = find(difTM >(dR/2));
bJumps = find(-difTM > (dR/2));

%nope - doesn't remove jump, just shortens it
%difTM(pRolInds)=difTM(pRolInds-1)+1;

difTM(fJumps)=(encVlt(fJumps)-dR)-encVlt(fJumps-1); %FROM ERIC

%nope - doesn't remove jump, just shortens it
%nRolInds=find(diff(encVlt)<(-dR/2));
%difTM(nRolInds)=difTM(nRolInds-1)-1;
difTM(bJumps)=(encVlt(bJumps)+dR)-encVlt(bJumps-1); %FROM ERIC

%changed for easier indexing
positionVec=smoothdata((-metPerInc)*cumsum(difTM),'lowess',200);

end



