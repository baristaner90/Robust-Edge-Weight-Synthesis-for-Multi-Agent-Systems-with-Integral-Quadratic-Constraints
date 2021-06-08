clear all
close all
clc
% % % Define systems
sys1        = ss([0,1;-2,-2],[0;1],[1,0],0);
sys2        = ss([0,1;-0.8,-2],[0;1],[1,0],0);

% Group of systems
Gco         = blkdiag(sys1,sys1,sys2);
nz1         = size(Gco.A,1)/size(sys1.A,1);
[A,B,C,D]   = ssdata(Gco);

% cooperative ready system
A_coop      = A;
B_coop      = [B,B,B];
C_coop      = [zeros(size(C));C;C];
D_coop      = zeros(size(B_coop,2),size(B_coop,2));
D_coop(1:nz1,2*nz1+1:3*nz1) = eye(nz1);
Gcoop       = ss(A_coop,B_coop,C_coop,D_coop); % {\bm S} in paper

% find sizes of the signals
deChn   =   size(Gco);
nz2     =   deChn(1);
nw2     =   deChn(2);
nz      =   nz2;
nw      =   nz;
nw1     =   size(Gcoop,2) - (nw2+nw);
nz1     =   size(Gcoop,1) - (nz2+nz);

sizes.nz1 = nz1;
sizes.nz2 = nz2;
sizes.nz  = nz;
sizes.nw1 = nw1;
sizes.nw2 = nw2;
sizes.nw  = nw;

alfa	= 0:0.1:1;
lam = [1 0];
% alfa	= 0.5;

for mm = 1:length(alfa)

Adj	= [0 , (alfa(mm)) , (1-alfa(mm));
        (alfa(mm)) , 0 , (1-alfa(mm));
        (1-alfa(mm)) , (alfa(mm)), 0];

%

        
[Hcoop, HElem]      =   hcoop(Gcoop,Adj,sizes,true);
gamma0_mat(mm)      =	norm(Hcoop,inf);
data0_mat(mm,:)     =   [alfa(mm),gamma0_mat(mm)];

mu      =   0.1;
gam_up	=	gamma0_mat(mm)*5;
gam_lw	=	0;
 

[solc, GMIN]	=	minGT_feasp_iqc(gam_up,gam_lw,Hcoop,HElem,sizes,lam);
gamma0_bt(mm)	=   solc.gamma;
data0_bt(mm,:)	=   [alfa(mm),gamma0_bt(mm)];

end

figure
plot(data0_bt(:,1),data0_bt(:,2));
% Phase 1: Find relative Interior for given Adjacency matrix
alf         =	0.3;
Adj0        =	[0    ,	alf	, 1-alf;
                 alf  ,   0	, 1-alf;
               1-alf  , alf	, 0    ];
           
[Hcoop0, HElem0]	=   hcoop(Gcoop,Adj0,sizes,true);

% Conventional method analysis elapsed cputime
lam0  = [1 0];
tim  = cputime;
% 
[sol0, GMIN]	=	minGT_feasp_iqc(gam_up,gam_lw,Hcoop0,HElem0,sizes,lam0);
% 
tim  =  cputime - tim;
% Prepare for second phase by combining Pi and creating equivalent HeElem

% Store data for second phase
solk0.Gk        =   sol0.gamma;
solk0.lamk      =   sol0.lam;
solk0.Xkvar     =   sol0.x;
solk0.alfa      =   alf;
solk0.Dak       =   Adj0; 
solk0.sys       =   Hcoop0;
solk0.sysElem   =	HElem0;
solk0.sizes     =   sizes;
solk0.Hek       =	sol0.He;
solk0.HeElemk   =	sol0.HeElem;
solk0.Psik      =	sol0.Psi;
solk0.PiElemk	=	sol0.PiElem;
solk0.Gcoop     =   Gcoop;
solk0.succeed   =   true;

% Second Phase: Augmented BMI for robust Adj Syn
solsit = solk0;
maxiter	= 10; 
for iter = 1:maxiter
        % Optimize Adjacency 
        [solsit] = seqlmiGlayer_iqc(solsit);
%         % Store data for second phase
%         solsit.Gk        =   solsseq.Gk;
%         solsit.lamk      =   solsseq.lamk;
%         solk0.Xkvar     =   solsseq.Xkvar;
%         solsit.alfa      =   solsseq.alphak;
%         solsit.Dak       =   solsseq.Dak; 
%         solsit.sys       =   solsit.sys;
%         solsit.sysElem   =	solsit.sysElem;
%         solsit.sizes     =   solsit.sizes;
%         solsit.Hek       =	solsit.Hek;
%         solsit.HeElemk   =	solsseq.HeElemk;
%         solsit.Psik      =	solsit.Psik;
%         solsit.PiElemk	=	solsseq.PiElemk;
%         solsit.Gcoop     =   solsseq.Gcoop;
%         solsit.succeed   =   true;
        solp2(iter)      =   solsit;
end
i_opt       =   length(solp2);
for i=1:i_opt
    l_opt(i,:)   =   solp2(i).lamk;
    g_opt(i,1)   =   solp2(i).Gk;
    a_str(i,1)   =   solp2(i).Dak(1,2);
end
% Optimize lambda
[Hcoopit, HElemit]	=   hcoop(Gcoop,solp2(end).Dak,sizes,true);
mu      =   0.1;
gam_up	=	solp2(end).Gk*2;
gam_lw	=	0;
[sols2, GMIN]	 =	minGT_feasp_iqc_p2(gam_up,gam_lw,Hcoopit,HElemit,sizes,solp2(end).lamk);

a_str = [a_str;solp2(end).alphak];
g_opt = [g_opt;sols2.gamma];
l_opt = [l_opt;sols2.lam];
%
alfa	= 0:0.1:1;
lam     = l_opt(end,:);
for mm = 1:length(alfa)

Adj	= [0 , (alfa(mm)) , (1-alfa(mm));
        (alfa(mm)) , 0 , (1-alfa(mm));
        (1-alfa(mm)) , (alfa(mm)), 0];

%

        
[Hcoop, HElem]      =   hcoop(Gcoop,Adj,sizes,true);
gamma0_mat(mm)      =	norm(Hcoop,inf);
data0_mat(mm,:)     =   [alfa(mm),gamma0_mat(mm)];

mu      =   0.1;
gam_up	=	gamma0_mat(mm)*5;
gam_lw	=	0;
 

[solc, GMIN]	=	minGT_feasp_iqc(gam_up,gam_lw,Hcoop,HElem,sizes,lam);
gamma_bt(mm)	=   solc.gamma;
data_bt(mm,:)	=   [alfa(mm),gamma_bt(mm)];
end

%
figure
% plot(data0_mat(:,1),data0_mat(:,2),'color','black','Marker','o','LineWidth',4);
% hold on
plot(data0_bt(:,1),data0_bt(:,2),'color','red','Marker','x','LineWidth',2,'MarkerSize',14);
hold on
plot(data_bt(:,1),data_bt(:,2),'color','blue','Marker','d','LineWidth',2,'MarkerSize',10);
plot(solk0.alfa, solk0.Gk,'color','Cyan','Marker','d','MarkerSize',14,'LineStyle','none','MarkerFaceColor','Magenta');
plot(a_str(1:end,1),g_opt(1:end,1),'color','black','Marker','s','MarkerSize',10,'LineStyle','--','LineWidth',1.1,'MarkerFaceColor','Green');
plot(a_str(end,1),g_opt(end,1),'color','black','Marker','d','MarkerSize',10,'LineStyle','--','LineWidth',1.1,'MarkerFaceColor','Red');
hold off
% title('Optimal Adjacency for given cooperative system with H_{\infty} Performance Criteria')
% Leg = legend( '||*||_{\infty} Norm along \alpha','\gamma w.r.t. \alpha for \lambda^{0}',...
%         '\gamma w.r.t. \alpha for \lambda^{opt}',...
%         '\gamma^{0}, \lambda=\lambda^{0}','\gamma^{t}','\gamma^{final}, \lambda=\lambda^{opt}');
Leg = legend('\gamma w.r.t. \alpha for \lambda^{0}',...
             '\gamma w.r.t. \alpha for \lambda^{opt}',...
             '\gamma^{0}, \lambda=\lambda^{0}','\gamma^{t}','\gamma^{final}, \lambda=\lambda^{opt}');
set(Leg, 'Location', 'Best');
xlabel('\alpha');
ylabel('\gamma');
% title('H_{\infty} Performance Plot')
% title('Synthesized \alpha vs \gamma')
xlim([data0_bt(1,1)-0.1 data0_bt(end,1)+0.1]);
ylim([0,max(data0_mat(:,2))+0.5])
grid on
savefig('robseqlmiiqc_het.fig')
saveas(gcf,'robseqlmiiqc_het','jpg')

%
[Hn, HnElem]	=   hcoop(Gcoop,Adj0,sizes,true);
An              =   ss(Adj0);
alfa_r          =   a_str(end);
Adjr            =	[0    ,	alfa_r      , 1-alfa_r;
                     alfa_r  ,   0      , 1-alfa_r;
                   1-alfa_r  , alfa_r	, 0    ];
[Hr, HrElem]	=   hcoop(Gcoop,Adjr,sizes,true);

% Run simulink file 
%%
T   = out.z.time;
zn  = out.z.Data(:,1:3);
zr  = out.z.Data(:,4:6);
figure
plot(T,zn(:,1),'Color','r','LineWidth',1.2);
hold on
plot(T,zr(:,1),'Color','k','LineWidth',1.2);
plot(T,zn(:,2),'Color','r','LineWidth',1.6,'LineStyle',':');
plot(T,zr(:,2),'Color','k','LineWidth',1.6,'LineStyle',':');
plot(T,zn(:,3),'Color','r','LineWidth',1.2,'LineStyle','-.');
plot(T,zr(:,3),'Color','k','LineWidth',1.2,'LineStyle','-.');
grid on
xlabel('Time, s')
ylabel('Output')
% title('Time Response Plot for Arbitrarily Injected Unit Impulse')
legend('Initial MAS G^{1}','Robust MAS G^{1}','Initial MAS G^{2}','Robust MAS G^{2}','Initial MAS G^{3}','Robust MAS G^{3}')