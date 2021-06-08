function [sols] = seqfeasCond_iqc(solk0,t)
sols.time  = cputime;
%  deploy parameters
sizes = solk0.sizes;
Gcoop = solk0.Gcoop;

%   Previous solutions
Gk      =	solk0.Gk;
lamk	=   solk0.lamk;
muk     =   sqrt(lamk);
Xkvar	=   solk0.Xkvar;
% Dak is not deployed since it is alrady in following system matrices
Aek     =   solk0.HeElemk.Ae;
Be1k    =   solk0.HeElemk.Be1;
Be2k    =   solk0.HeElemk.Be2;
Zp1k    =   solk0.HeElemk.Zp1;
Zp2k    =   solk0.HeElemk.Zp2;
Zek     =   solk0.HeElemk.Z2;
% Sizes 
nx      =   size(Aek,1);
% np      =   size(Zp1k,1);
nw1     =   sizes.nw1;
nw2     =   sizes.nw2;
npi     =   solk0.PiElemk.npi;
%   Correlation Matrix, Pp=vp'*Qp*vp
%   vp*(vp'*Qp*vp)*vp' = [0, I; I, 0]
Qp      =   [zeros(nw1), eye(nw1);
             eye(nw1), zeros(nw1)];
[vp,ep] =   eig(Qp);
vp      =	[vp(:,nw1+1:end),vp(:,1:nw1)];
ep      =   blkdiag(eye(nw1),-eye(nw1));
% Write modified Zp1k and Zp2k
Zp1kk    =   vp*Zp1k*muk(1);
Zp2kk    =   vp*Zp2k*muk(2);

% Partition Zp1k and Zp2k
% Zp1k  =   [Zp11k;Zp12k]
% Zp2k  =   [Zp21k;Zp22k]
Zp11k   =   Zp1kk(1:nw1,:);
Zp12k   =   Zp1kk(nw1+1:end,:);
Zp21k   =   Zp2kk(1:nw1,:);
Zp22k   =   Zp2kk(nw1+1:end,:);

% Z1k = Zp11k'*Zp12k and Z2k = Zp21k'*Zp22k
Rk     =   [Zp11k; Zp21k];
Lk     =   [Zp12k; Zp22k];

% Create Pk Qk using vp
Qk  =   zeros(size(Zek));  Qk(:,end-nw2+1:end)=eye(nw2);
Pk  =   (1/Gk)*Zek;

% Transform Pk Qk using vp
PQk =   [Pk;Qk];
PQkk=   vp*PQk;
% Seperate \bar{Pk} and \bar{Qk}
Pbk =   PQkk(1:nw2,:);
Qbk =   PQkk(nw2+1:end,:);

Xk	=	[Xkvar, zeros(size([Be1k, Be2k]))];
Yk  =   [Aek, Be1k, Be2k];

Mk  =   [Xk;Pbk;Rk;Yk;Qbk;Lk];
%   Linearization takes place is the following code.

%   Variables of the problem
Xvar	=   sdpvar(nx, nx, 'symmetric');
da1     =   sdpvar(1,1);
Da      =   [0,da1,1-da1;
             da1,0,1-da1;
             1-da1,da1,0]; 

% Extended system is created with variables
[~, HElem]      =   hcoop(Gcoop,Da,sizes,false);
[Psi,PiElem,~]	=	creaiqc(sizes);
[~,HeElemvar]	=	hext([],HElem,Psi,PiElem,sizes,false);
% 
Avar    =   HeElemvar.Ae;
Bvar    =   [HeElemvar.Be1,HeElemvar.Be2];
Zevar   =   HeElemvar.Z2;
Zp1var  =   HeElemvar.Zp1;
Zp2var  =   HeElemvar.Zp2;
Zp1var	=   vp*Zp1var*muk(1);
Zp2var	=   vp*Zp2var*muk(2);

% Partition
Zp11var   =   Zp1var(1:nw1,:);
Zp12var   =   Zp1var(nw1+1:end,:);
Zp21var   =   Zp2var(1:nw1,:);
Zp22var   =   Zp2var(nw1+1:end,:);

% Compose variable R and L 
R       =   [Zp11var; Zp21var];
L       =   [Zp12var; Zp22var];
% Create P and Q
Q       =   zeros(size(Zek));  Q(:,end-nw2+1:end)=eye(nw2);
P       =   (1/Gk)*Zevar;
% Transform P Q using vp
PQ   =   [P;Q];
PQvv =   vp*PQk;
% Seperate \bar{Pk} and \bar{Qk}
Pb =   PQvv(1:nw2,:);
Qb =   PQvv(nw2+1:end,:);

X       =	[Xvar, zeros(size(Bvar))];
Y       =   [Avar, Bvar];

Mb  =   [X;Pb;R;Y;Qb;L];

% Make Qb to collect everything together
Ib      =   eye(size(Mb,1)/2);
Ob      =   zeros(size(Ib));
Qb      =   [Ob, Ib; Ib, Ob];
[vb,~]	=   eig(Qb);
% Qb1     =   vb*[Ob; Ib]*[Ob; Ib]'*vb';
Qb2     =   vb*[Ib; Ob]*[Ib; Ob]'*vb';

LM12    =   Mb'*Qb2*Mk + Mk'*Qb2*Mb - Mk'*Qb2*Mk;

Phi11	=   -Ib;
Phi12	=   ([Ob; Ib]')*(vb')*Mb;
Phi21	=   Phi12';
Phi22	=   - LM12 - t*eye(size(LM12));

Phi     =   [Phi11, Phi12; Phi21, Phi22];

F       =   [Phi<=0; Xvar>=0; 1>=da1>=0];
% F       =   [Phi<=0; Xvar>=0; 1>=da1>=0; Ix<=Xvar<=gc*Ix; lam(1)>=0; lam(2)>=0];
opt     =   sdpsettings('verbose',0,'warning',0);
diag    =   optimize(F,[],opt);
%   Store Solutions
sols.succeed     =   false;
sols.tk          =   [];
sols.Xkvar       =   [];
sols.Dak         =	 [];
sols.Gk          =   [];
sols.alphak      =   [];
sols.MaxEigenLMI =   [];
sols.MinEigenLMI =   [];
sols.MaxEigenX   =   [];
sols.MinEigenX   =   [];
sols.FrobRadLim	 =   [];
sols.FrobRadAct  =   [];
sols.FrobRadSat  =   [];
sols.sizes       =   [];
sols.Gcoop       =   [];
sols.lamk        =   [];
sols.HeElemk     =	 [];
sols.PiElemk     =   [];
if diag.problem~=0
    sols.time        = cputime - sols.time;
elseif diag.problem==0
    if ~isempty(value(Xvar)) && ~isempty(value(Da))
        if t<0
            xval    =   value(Xvar);
            Daval   =   value(Da);
            lamval  =   value(lamk);
            % Extended system is created with variables
            [~, HElemval]	=   hcoop(Gcoop,Daval,sizes,false);
            [~, Heval]      =	hext([],HElemval,Psi,PiElem,sizes,false);
            Aval    =   Heval.Ae;
            Bval    =   [Heval.Be1,Heval.Be2];
            Zeval   =   Heval.Z2;
            Zp1val  =   Heval.Zp1;
            Zp2val  =   Heval.Zp2;
            Xval	=	[xval, zeros(size(Bval))];
            Yval	=   [Aval, Bval];
            Zval	=   (1/(Gk^2))*(Zeval')*Zeval + blkdiag(zeros(nx+nw1),-eye(nw2));
            lmicheck	=   Xval'*Yval+Yval'*Xval+Zval+lamval(1)*Zp1val'*ep*Zp1val+lamval(2)*Zp2val'*Qp*Zp2val;
            max_eigenLMI	=   max(eig(lmicheck));
            min_eigenLMI	=   min(eig(lmicheck));
            if max_eigenLMI < 0 && min_eigenLMI < 0
                max_eigenX  =   max(eig(xval));
                min_eigenX  =   min(eig(xval));
                EigX        =	eig(xval);
                EucNorm     =   norm(EigX,2);
%                 FroNorm     =   norm(xval,'fro');
                sols.succeed     = true;
                sols.tk          = t;
                sols.Xkvar       = xval;
                sols.Gk          = Gk;
                sols.Dak         = Daval;
                sols.alphak      = value(da1);
                sols.sizes       = sizes;
                sols.Gcoop       = Gcoop;
                sols.lamk        = lamval;
                sols.HeElemk     = Heval;
                sols.PiElemk     = PiElem; 
                sols.MaxEigenLMI = max_eigenLMI;
                sols.MinEigenLMI = min_eigenLMI;
                sols.MaxEigenX   = max_eigenX;
                sols.MinEigenX   = min_eigenX;
                sols.FrobRadLim  = inf;
                sols.FrobRadAct  = EucNorm;
                sols.FrobRadSat  = (EucNorm/sols.FrobRadLim)*100;
            end
        end
    end
	sols.time        = cputime - sols.time;       
end
end