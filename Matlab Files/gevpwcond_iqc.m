function [sol] = gevpwcond_iqc(g,fr,He,HeElem,lam)
% Deploy HeElem
Ae      =   HeElem.Ae;
Be1     =   HeElem.Be1;
Be2     =   HeElem.Be2;
Ce1     =   HeElem.Ce1;
Ce2     =   HeElem.Ce2;
De11    =   HeElem.De11;
De12    =   HeElem.De12;
De21    =   HeElem.De21;
De22    =   HeElem.De22;
Z2      =   HeElem.Z2;
Zp      =   HeElem.Zp;
Zp1     =   HeElem.Zp1;
Zp2     =   HeElem.Zp2;
nx      =   size(Ae,1);
nw1     =   size(Be1,2);
nw2     =   size(Be2,2);
G       =   g;
G2inv   =   (G^2)^(-1);
IG      =   blkdiag(zeros(nx+nw1),-eye(nw2));
Ie      =   [eye(nx);zeros(size([Be1, Be2]'))];
lam1    =   lam(1);
lam2    =   lam(2);
Ppsi    =   blkdiag(eye(nw2),-eye(nw2));

setlmis([]);
[X,nX,sX] = lmivar(1, [nx,1]);
% for i=1:2
%     [lam(i),nlam(i),slam(i)] = lmivar(1, [1,1]);
% end

nlmi	= 1;
lmiterm([nlmi 1 1 X], Ie, [Ae,Be1,Be2], 's');
lmiterm([nlmi 1 1 0], G2inv*Z2'*Z2);
lmiterm([nlmi 1 1 0], IG);
lmiterm([nlmi 1 1 0], lam1*Zp1'*Ppsi*Zp1);
% lmiterm([nlmi 1 1 lam(1)], 1,Zp1'*Ppsi*Zp1);
lmiterm([nlmi 1 1 0], lam2*Zp2'*Ppsi*Zp2);
% lmiterm([nlmi 1 1 lam(2)], 1,Zp2'*Ppsi*Zp2);

nlmi = nlmi + 1;
lmiterm([-nlmi 1 1 X], 1, 1);
% nlmi = nlmi + 1;
% lmiterm([-nlmi 1 1 lam(1)], 1, 1);
% nlmi = nlmi + 1;
% lmiterm([-nlmi 1 1 lam(2)], 1, 1);

lmis = getlmis;
options = [0,0,fr,0,1];
% clc
[tmin,xfeas] = feasp(lmis,options);
sol.succeed     = false;
sol.t           = [];
sol.x           = [];
sol.lam         = [];
sol.gamma       = [];
sol.MaxEigenLMI = [];
sol.MinEigenLMI = [];
sol.MaxEigenX   = [];
sol.MinEigenX   = [];
sol.FrobRadLim	= [];
sol.FrobRadAct  = [];
sol.FrobRadSat  = [];
if ~isempty(xfeas)
    if tmin<0
        xval        =   dec2mat(lmis, xfeas, X);
%         lamval1      =   dec2mat(lmis, xfeas, lam(1));
%         lamval2      =   dec2mat(lmis, xfeas, lam(2));
%         lamval      =   [lamval1, lamval2];
          lamval    =   [lam1, lam2];
        lmicheck    =   Ie*xval*[Ae,Be1,Be2]+(Ie*xval*[Ae,Be1,Be2])' + ...
                        G2inv*(Z2')*Z2 +IG + lamval(1)*(Zp1')*Ppsi*Zp1 + ...
                        lamval(2)*(Zp2')*Ppsi*Zp2;
%         lmicheck    =   Ie*xval*[Ae,Be1,Be2]+(Ie*xval*[Ae,Be1,Be2])' + ...
%                         G2inv*(Z2')*Z2 +IG + lamval1*(Zp1')*Ppsi*Zp1 + ...
%                         lamval2*(Zp2')*Ppsi*Zp2;
        max_eigenLMI  =   max(eig(lmicheck));
        min_eigenLMI  =   min(eig(lmicheck));
        if max_eigenLMI < 0 && min_eigenLMI < 0
            max_eigenX  =   max(eig(xval));
            min_eigenX  =   min(eig(xval));
            EigX        =	eig(xval);
            EucNorm     =   norm(EigX,2);
%             FroNorm     =   norm(xval,'fro');
            sol.succeed     = true;
            sol.t           = tmin;
            sol.x           = xval;
            sol.lam         = lamval; %[lam1 lam2];
            sol.gamma       = G;
            sol.MaxEigenLMI = max_eigenLMI;
            sol.MinEigenLMI = min_eigenLMI;
            sol.MaxEigenX   = max_eigenX;
            sol.MinEigenX   = min_eigenX;
            sol.FrobRadLim	= fr;
            sol.FrobRadAct  = EucNorm;
            sol.FrobRadSat  = (EucNorm/fr)*100;
        end
    else
    end
end
    
end