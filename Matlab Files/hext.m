function [He,HeElem] = hext(Hcoop,HElem,Psi,PiElem,sizes,flag)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Deploy sizes
nz1 = sizes.nz1;
nz2 = sizes.nz2;
nz  = sizes.nz;
nw1 = sizes.nw1;
nw2 = sizes.nw2;
nw  = sizes.nw;
%  Deploy HElem
A       = HElem.A;
B1      = HElem.B1;
B2      = HElem.B2;
C1      = HElem.C1;
C2      = HElem.C2;
D11     = HElem.D11;
D12     = HElem.D12;
D21     = HElem.D21;
D22     = HElem.D22;
%  Deploy PiElem
ap      = PiElem.ap;
bpz1	= PiElem.bpz1;
bpw1	= PiElem.bpw1;
cp      = PiElem.cp;
dpz1	= PiElem.dpz1;
dpw1	= PiElem.dpw1;
npi     = PiElem.npi;
% 
Ae      =   [A, zeros(size(A,1),size(ap,2));bpz1*C1, ap];
Be1     =   [B1;bpz1*D11+bpw1];
Be2     =   [B2;bpz1*D12];
Ce1     =   [bpz1*C1, cp];
if isempty(Ce1)
    Ce1 = zeros(nw1+nz1,size(A,2));
end
Ce2     =   [C2, zeros(size(C2,1),size(cp,2))];
De11    =   [dpz1*D11+dpw1];
De12    =   [dpz1*D12];
De21    =   [D21];
De22    =   [D22];
if flag
    He      =   ss(Ae,[Be1,Be2],[Ce1;Ce2],[De11, De12;De21, De22]);
else
    He      =   [];
end
HeElem.Ae   =   Ae;
HeElem.Be1  =   Be1;
HeElem.Be2  =   Be2;
HeElem.Ce1  =   Ce1;
HeElem.Ce2  =   Ce2;
HeElem.De11 =   De11;
HeElem.De12 =   De12;
HeElem.De21 =   De21;
HeElem.De22 =   De22;
HeElem.Z2   =   [Ce2,De21,De22];
HeElem.Zp   =   [Ce1,De11,De12];
HeElem.Zp1  =   HeElem.Zp(1:+nz1+nw1,:);
HeElem.Zp2  =   HeElem.Zp(1+nz1+nw1:end,:);
end

