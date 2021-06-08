function [sols] = seqlmiGlayer_iqc(solk0)
%   Options For gamma minimization
gam_lw	=   0;
gam_up	=   solk0.Gk*2;
gam_tol	=   0.01;
gam_err	=   inf;
gam_old	=   0;
max_iter=   20;

gam_cnt	=   1;
while true
    yalmip('clear')
    gam_try	=   (gam_lw + gam_up)*0.5
    solk0.Gk = gam_try;
    ttrial = -0.00001;
    [sols]	=	seqfeasCond_iqc(solk0,ttrial);
    if sols.succeed
        gam_up	=	gam_try;
    else
        gam_lw	=	gam_try;
        [sols]	=	solk0;
        sols.succeed    =   false;
    end
    gam_err   =   abs(gam_old - gam_try);
    if sols.succeed && gam_err<gam_tol
        sols.gam_err      =   gam_err;
        sols.gam_iter     =   gam_cnt;
        sols.Gk           =   gam_try;
        break
    elseif gam_cnt>max_iter
        break
    end
    gam_cnt           =   gam_cnt +  1;
    gam_old           =   gam_try;    
end
end