%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Flavio Luiz Cardoso-Ribeiro: http://github.com/flavioluiz/     %
%           ISAE-Supaero   / Instituto Tecnologico de Aeronautica  %
%                               CNPq   - Brazil                    %
%                                                                  %
%    This project is part of ANR Project HAMECMOPSYS:              %
%                  https://hamecmopsys.ens2m.fr/                   %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fluid = rotsaintvenant(numelements,FILLING)
    % semi-discretization saint-venant under rotation
    % Xp = J Q X + Bu
    % y = B' QX + Du
    %
    % X = [alp1;alp2;pf;theta]
    %
    % J: (2*n+2) x (2*n+2)
    % Q: (2*n+2) x (2*n+2)
    % 
    %      u = [e1(L)  e2(0) Mext]        %                                 
    %      y = [-e2(L) e1(0) thetap]
    %  u = [FL,v0, Mext]   
    %  y = [-vL,F0, thetap]
    n = numelements;
    [pslosh p prigid] = dataexperiment(FILLING); % look dataexperiment.m file:
    rho = pslosh.rho;
    g = pslosh.g;
    L = pslosh.a;
    b = pslosh.b;
    hbar = pslosh.h;
    e = ones(n,1);
    diags = zeros(n,n+1);
    alp = 1;
    for i = 1:n+1
        k = (n+2-i);
        if k == 1 
            diags(:,i) = -e*(1/alp);
        else
            diags(:,i) = e*(1/alp)^2*((alp-1)/alp)^(k-2);
        end    
    end
    M = spdiags(diags, -n:0, n, n);
    J = [M*0 M; -transpose(M) M*0];

    Bq = zeros(n,1);
    Bp = zeros(n,1);
    for i = 1:n
        Bq(i) = 1/alp*((alp-1)/alp)^(i-1);
        k = n+1-i;
        Bp(i) = -1/alp*((alp-1)/alp)^(k-1);
    end
    dx = L/(n);
    B = [0*Bq Bq; Bp 0*Bp];
    B = [B zeros(2*n,1)];
    B = [B;...
         0,0,1;...
         0,0,0];
    
    D = [0, ((alp-1)/alp)^n;-((alp-1)/alp)^n 0];
    D = blkdiag(D,0);
    % new J matrix with rigid-body part:
    Jrb = [0,-1;1,0];
    J = blkdiag(J,Jrb);
    
    % defining Q for saint-venant equation (no rotation)    
    onesd = eye(n);
    Q1 = onesd*rho*b*g*hbar^2/dx; %only sloshing part
    Q2 = onesd*1/(rho*hbar*b)/dx;
    Q = blkdiag(Q1,Q2);
        
    Q1theta = zeros(n,1);
    bz = @(i)(i*dx-L/2);
    for i = 1:n
        Q1theta(i,1) = rho*b*g*hbar*(bz(i)^2-bz(i-1)^2)/2/dx;%Q1theta(i,1) = rho*b*g*hbar*(2*bz(i)+dx);
    end
    eig(Q)
    Q = [Q,zeros(2*n,1),[Q1theta;Q1theta*0]];
    Q = [Q; zeros(1,n*2+2); Q1theta',Q1theta'*0,0,0];        
    
    If = prigid.Ifluid;
    Q(end-1,end-1) = 1/If;
    Q(end,end) = 0;
    eig(Q)
%     size(J)
%     size(Q)
%     val = sort(damp(eig(J*Q)));
%     val(1:5)
%     pi*sqrt(g*hbar)/L*[1:2:7]/2

    fluid.J = J;
    fluid.Q = Q;
    fluid.B = B;
    fluid.D = D;
%     M = [1,0,0;
%          0,0,0;
%          0,0,1];
%     N = [0,0,0;
%          0,1,0;
%          0,0,0];
%     A = [J*Q B; M*transpose(B)*Q M*D+N];
%     E = eye(size(A));
%     E(end-2:end,end-2:end) = 0;
%     sort(damp(eig(A,E))/2/pi)
end