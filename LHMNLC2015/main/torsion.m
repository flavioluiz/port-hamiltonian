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

function torsionbeam = torsion(numelem)
    % semi-discretization of wave equation: torsion beam
    % Xd = JQ X + B u
    % y = B' QX + D u
    %
    % u = [e1(L)  e2(0)] [  -ML   -thetad0] (torsion moment at free tip)
    % y = [-e2(L)  e1(0)] [thetadL M0] (angle of torsion rate at free tip)
    
    n = numelem;

    e = ones(n,1);
    diags = zeros(n,n+1);
    alp = 0.5;
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

    B = [0*Bq Bq; Bp 0*Bp];

    D = -[0, ((alp-1)/alp)^n;-((alp-1)/alp)^n 0];

    
    [dum beamdata dum2] = dataexperiment(0); %load beam data
    % defining Q for torsion equation:
    E = beamdata.E; %Young modulus
    mu = beamdata.mu;  %Poisson ratio
    rho = beamdata.rho; %density
    L = beamdata.L; %length
    h = beamdata.h;  %thickness       +- 0.00005
    l = beamdata.l;   %width           +- 0.001
    I = l*h^3/12;
    I2 = l^3*h/12;    
    JJ = 0.33*l*h^3;
    G = E/2/(1+mu);
    Ip = (I+I2)*rho;    
    GJ = G*JJ;
    
    dx = L/(n);
    onesd = eye(n);
    Q = blkdiag(onesd/dx*GJ,onesd*1/(Ip)/dx);
    %val = sort(damp(eig(J*Q)));
    %val(1:2:9)'
    %sqrt(GJ/Ip)/L*[1:2:7]/2  /2
    torsionbeam.J = J;
    torsionbeam.Q = Q;
    torsionbeam.B = B(:,1);
    torsionbeam.D = D(1,1);
    
%     [eigvec eigval] = eig(J*Q);
%     M = [1,0;0,1];
%     N = [0,0;0,0];
%     A = [J*Q, B; M*transpose(B)*Q M*D+N];
%     E = eye(size(A));
%     E(end-1:end,end-1:end) = 0;
%     autoval = eig(A,E);
%     an = sort(damp(autoval)/2/pi);
%     an(1:2:6)
%     [autoval,sor] = sort(damp(diag(eigval))/2/pi);
%     autoval(1:2:6)
end