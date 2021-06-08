function [solc, GMIN] = minGT_feasp_iqc(gam_up,gam_lw,Hcoop,HElem,sizes,lam)
%   Options For gamma minimization

nz1     =   sizes.nz1;
gam_tol	=   0.001;
gam_err	=   inf;
gam_old	=   0;

gam_cnt	=   1;
while true
    
    gam_try     =   (gam_lw + gam_up)*0.5;
    % Define IQCs
    Psis        = pcon_2(nz1);
%     [Pi_lam,Psi,P_lam,PiElem] = combPi(Psis, lam, sizes);
    [Psi,PiElem,~]      =	creaiqc(sizes);
    [He,HeElem]         =	hext(Hcoop,HElem,Psi,PiElem,sizes,true);
    % [solc,~]      =	minTLayer(gam_try,Hcl_inf);
    [solc]              =   gevpCond_iqc(gam_try,He,HeElem,lam);
    if solc.succeed
        gam_up      =   gam_try
        gam_str(gam_cnt)  =   gam_try;   
    else
        gam_lw	=  gam_try;
    end
    gam_err	=   abs(gam_old - gam_try);
    if gam_err<gam_tol && solc.succeed
        GMIN.gam_err	=   gam_err;
        GMIN.gam_iter	=   gam_cnt;
        GMIN.gam        =   gam_try;
        break
    else
        gam_cnt         =   gam_cnt +  1;
        gam_old         =   gam_try;
    end
end
i = length(gam_str);
while true
    if gam_str(i)~=0
        GAM     =   gam_str(i);
        break
    else
        i=i-1;
    end
end

[solc]      = gevpCond_iqc(GAM,He,HeElem,lam);
solc.He     = He;
solc.HeElem = HeElem;
solc.Psi    = Psi;
solc.PiElem	= PiElem;

fprintf("==============================")
fprintf("\n")
fprintf("End of OPT")
fprintf("\n")
fprintf("==============================")
fprintf("\n")
end