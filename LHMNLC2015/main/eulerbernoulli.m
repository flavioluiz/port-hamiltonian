%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%																   %
%   Flavio Luiz Cardoso-Ribeiro: https://github.com/flavioluiz/    %
%			ISAE-Supaero   / Instituto Tecnologico de Aeronautica  %
%								CNPq   - Brazil 				   %
%																   %
%    This project is part of ANR Project HAMECMOPSYS:              %
%					https://hamecmopsys.ens2m.fr/ 				   %
%																   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ebbeam = eulerbernoulli(numelem)
% semi-discretization of second-order equation (like Euler-Bernoulli):
% Xp = JQX + Bu
% y = B'*QX + Du
% x = [];
% where 
%       u = [-e1^B(L) -e1d^B(L) -e2(0) -e2d(0)] (-moment, force in the "free" tip)
%       y = [-e2d^B(L) e2^B(L)  -e1d(0) e1(0)]  (angle rate, speed in the tip)

%       u = [mL fL wd0 thetad0] (-moment, force in the "free" tip)
%       y = [-thetadL wdL -f0 m0]  (angle rate, speed in the tip)
%
    [dum beamdata dum2] = dataexperiment(0);
    L = beamdata.L;

    EI = beamdata.E*beamdata.l*(beamdata.h^3)/12;
    rho = beamdata.rho*beamdata.l*beamdata.h;

    n = numelem;
    e = ones(n,1);
    diags = zeros(n,n+1);
    alp = 0.5;
    
    dx = L/(n);
    alpb = 1/dx;
    for i = 1:n+1
        k = (n+2-i);
        if k == 1 
            diags(:,i) = -e*(alpb/alp^2);
        else
            diags(:,i) = e*alpb*(2*alp-k)*(alp-1)^(k-3)/alp^(k+1);
        end    
    end
    M = spdiags(diags, -n:0, n, n);
    J = [M*0 M; -transpose(M) M*0];

    B1 = zeros(n,2);
    B2 = zeros(n,2);
    for i = 1:n
        B1(i,:) = [alpb*(alp-i)*(alp-1)^(i-2)/alp^(i+1), (alp-1)^(i-1)/alp^i ];
        k = n+1-i;
        B2(i,:) = [-alpb*(alp-k)*(alp-1)^(k-2)/alp^(k+1), (alp-1)^(k-1)/alp^k ];
    end

    B = [B1*0 B1; B2 B2*0];

       %check the sign of D!
    D = [0, 0, -n * alpb*(alp-1)^(n-1)/(alp)^(n+1), -((alp-1)/alp)^n;
         0, 0,  ((alp-1)/alp)^n, 0;
         n * alpb*(alp-1)^(n-1)/alp^(n+1), -((alp-1)/alp)^n, 0, 0;
         ((alp-1)/alp)^n, 0, 0,0];

    % defining Q for Euler-Bernoulli equation:

    
    onesd = eye(n);
    Q = blkdiag(onesd*EI/dx,onesd*1/(rho)/dx);

    %val = sort(damp(eig(J*Q)));
    %val(1:2:10)/2/pi;
    %pi*sqrt(g*hbar)/L*[1:2:7]/2

    ebbeam.J = J;
    ebbeam.Q = Q;
    ebbeam.B = B(:,1:2);
    ebbeam.D = D(1:2,1:2);    
end