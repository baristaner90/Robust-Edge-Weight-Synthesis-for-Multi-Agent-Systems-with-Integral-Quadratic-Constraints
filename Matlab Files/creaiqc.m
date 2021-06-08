function [Psi,Elem,Psis] = creaiqc(sizes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Deploy sizes
nz1 = sizes.nz1;
nw1 = sizes.nw1;


Psis        = pcon_2(nz1);
% Evaluate number of IQCs
npi = size(Psis,2);
% Construct the concetanated iqc system having z1 and w1 as inputs
Psi_con	=   [];
Psi_in	=   blkdiag(eye(nz1),eye(nw1));
for i=1:npi
    Psi_con	=   append(Psi_con,Psis{i});
    if i>1
        Psi_in	=   cat(1,Psi_in, eye(nz1+nw1));
    end
end
Psi_in  =   ss([],[],[],Psi_in);
Psi     =   series(Psi_in,Psi_con);
[ap,bp,cp,dp]   =   ssdata(Psi);
bpz1    =   bp(:,1:nz1);
bpw1    =   bp(:,nz1+1:end);
dpz1    =   dp(:,1:nz1);
dpw1    =   dp(:,nz1+1:end);
Elem.ap     =   ap;
Elem.bpz1   =   bpz1;
Elem.bpw1   =   bpw1;
Elem.cp     =   cp;
Elem.dpz1   =   dpz1;
Elem.dpw1   =   dpw1;
Elem.npi    =   npi;

end

