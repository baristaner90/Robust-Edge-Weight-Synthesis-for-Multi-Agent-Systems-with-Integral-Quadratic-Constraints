function [sol] = gevpCond_iqc_p2(g,He,HeElem,lam)
%	radius for conditioning is defined in feasp by 3rd option entry which
%	dictates the upper bound for Euclidean Norm of decision variable
fr	= 0;                 
[sol] = gevpwcond_iqc_p2(g,fr,He,HeElem,lam);
%   after having solution for unbounded case for the norm we can set that
%   value as an upper bound for the fr. Lower bound will still be 0.
if sol.succeed
    fr_lw   = fr;
    fr_up   = sol.FrobRadAct;
    tol         =   0.0001;
    fr_err      =   inf;
    fr_old      =   0;
    cnt         =   1;
    while  fr_err > tol
        fr_try   =  (fr_lw + fr_up)*0.5;
        [sol] = gevpwcond_iqc_p2(g,fr_try,He,HeElem,lam);
        if sol.succeed
            fr_up   =   fr_try;
        else
            fr_lw   =   fr_try;
        end
        fr_err  =   abs(fr_old - fr_try);
        if fr_err < tol
            break
        else
            cnt     =   cnt + 1;
            fr_old  =   fr_try;
        end
    end
    fr_margin   =   1;
    FR  =   fr_try+fr_margin;
    [sol] = gevpwcond_iqc_p2(g,FR,He,HeElem,lam);
    sol.FR  =   FR;
else
    fprintf("==============================")
    fprintf("\n")
	fprintf("NO SOLUTION at FR level")
    fprintf("\n")
    fprintf("==============================")
    fprintf("\n")

end
end