function [Pcon] = pcon_2(N)
% provide uncertain dynamics for multiagent system with interconnection topology
% as a constraint using IQCs (there will be 2 IQS)
Psi1 = ss(eye(2));

Hsat = ss(tf(1,[1 1]));
Pisat2 = [0, (1 + Hsat); (1 + Hsat'), -(2 + Hsat + Hsat')];
Psisat2 = jspec(Pisat2);

Pisat3 = [0 1; 1 -2];
Psisat3 = jspec(Pisat3);

rho = 1;
ub = 1;
Hslow = ss(tf([1 2],[1 1]));
Pisat4 = [(1+1)*(Hslow'*Hslow+1), 0; 0, -Hslow'*Hslow];
Psisat4 = jspec(Pisat4);

w=1.5; Theta0 = 1;
Psi0 = (w^2+0.08*w^4)/(1+0.13*w^2+0.02*w^4);
Pisat5 = [Psi0, 0;0, -1];
Psisat5 = jspec(Pisat5);

Psi1_11 = Psi1(1,1);
Psi1_12 = Psi1(1,2);
Psi1_21 = Psi1(2,1);
Psi1_22 = Psi1(2,2);

Psi2_11 = Psisat2(1,1);
Psi2_12 = Psisat2(1,2);
Psi2_21 = Psisat2(2,1);
Psi2_22 = Psisat2(2,2);

Psi3_11 = Psisat3(1,1);
Psi3_12 = Psisat3(1,2);
Psi3_21 = Psisat3(2,1);
Psi3_22 = Psisat3(2,2);

Psi4_11 = Psisat4(1,1);
Psi4_12 = Psisat4(1,2);
Psi4_21 = Psisat4(2,1);
Psi4_22 = Psisat4(2,2);

Psi5_11 = Psisat5(1,1);
Psi5_12 = Psisat5(1,2);
Psi5_21 = Psisat5(2,1);
Psi5_22 = Psisat5(2,2);

for i=1:N-1
    Psi1_11 = blkdiag(Psi1_11,Psi1(1,1));
    Psi1_12 = blkdiag(Psi1_12,Psi1(1,2));
    Psi1_21 = blkdiag(Psi1_21,Psi1(2,1));
    Psi1_22 = blkdiag(Psi1_22,Psi1(2,2));
    
    Psi2_11 = blkdiag(Psi2_11,Psisat2(1,1));
    Psi2_12 = blkdiag(Psi2_12,Psisat2(1,2));
    Psi2_21 = blkdiag(Psi2_21,Psisat2(2,1));
    Psi2_22 = blkdiag(Psi2_22,Psisat2(2,2));
    
    Psi3_11 = blkdiag(Psi3_11,Psisat3(1,1));
    Psi3_12 = blkdiag(Psi3_12,Psisat3(1,2));
    Psi3_21 = blkdiag(Psi3_21,Psisat3(2,1));
    Psi3_22 = blkdiag(Psi3_22,Psisat3(2,2));

    Psi4_11 = blkdiag(Psi4_11,Psisat4(1,1));
    Psi4_12 = blkdiag(Psi4_12,Psisat4(1,2));
    Psi4_21 = blkdiag(Psi4_21,Psisat4(2,1));
    Psi4_22 = blkdiag(Psi4_22,Psisat4(2,2));    
    
    Psi5_11 = blkdiag(Psi5_11,Psisat5(1,1));
    Psi5_12 = blkdiag(Psi5_12,Psisat5(1,2));
    Psi5_21 = blkdiag(Psi5_21,Psisat5(2,1));
    Psi5_22 = blkdiag(Psi5_22,Psisat5(2,2));   
end
Pcon{1} = [Psi1_11, Psi1_12; Psi1_21, Psi1_22];
Pcon{2} = [Psi2_11, Psi2_12; Psi2_21, Psi2_22];
% Pcon{3} = [Psi3_11, Psi3_12; Psi3_21, Psi3_22];
% Pcon{4} = [Psi4_11, Psi4_12; Psi4_21, Psi4_22];
% Pcon{5} = [Psi5_11, Psi5_12; Psi5_21, Psi5_22];
end