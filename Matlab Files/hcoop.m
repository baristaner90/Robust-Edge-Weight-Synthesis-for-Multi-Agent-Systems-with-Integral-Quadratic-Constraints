function [Hcoop, Elem] = hcoop(Gcoop,Da,sizes,flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Deploy sizes
nz1 = sizes.nz1;
nz2 = sizes.nz2;
nz  = sizes.nz;
nw1 = sizes.nw1;
nw2 = sizes.nw2;
nw  = sizes.nw;
% Make proper fractioning of the system
[a,b,c,d]   =   ssdata(Gcoop);
b1  =   b(:,1:nw1);
b2  =   b(:,nw1+1:nw1+nw2);
b   =   b(:,nw1+nw2+1:end);
c1  =   c(1:nz1,:);
c2  =   c(nz1+1:nz1+nz2,:);
c   =   c(nz1+nz2+1:end,:);
d11	=   d(1:nz1,1:nw1);
d12 =   d(1:nz1,nw1+1:nw1+nw2);
e1  =   d(1:nz1,nw1+nw2+1:end);
d21 =   d(nz1+1:nz1+nz2,1:nw1);
d22 =   d(nz1+1:nz1+nz2,nw1+1:nw1+nw2);
e2  =   d(nz1+1:nz1+nz2,nw1+nw2+1:end);
f1  =   d(nz1+nz2+1:end,1:nw1);
f2  =   d(nz1+nz2+1:end,nw1+1:nw1+nw2);
e   =   d(nz1+nz2+1:end,nw1+nw2+1:end);

A   =   [a+b*Da*c];
B1  =   [b1+b*Da*f1];
B2  =   [b2+b*Da*f2];
C1  =   [c1+e1*Da*c];
C2  =   [c2+e2*Da*c];
D11 =   [d11+e1*Da*f1];
D12 =   [d12+e1*Da*f2];
D21 =   [d21+e2*Da*f1];
D22 =   [d22+e2*Da*f2];
if flag
    Hcoop       = ss(A,[B1,B2],[C1;C2],[D11,D12;D21,D22]);
else
    Hcoop = [];
end

Elem.A      = A;
Elem.B1     = B1;
Elem.B2     = B2;
Elem.C1     = C1;
Elem.C2     = C2;
Elem.D11    = D11;
Elem.D12    = D12;
Elem.D21    = D21;
Elem.D22    = D22;
end

