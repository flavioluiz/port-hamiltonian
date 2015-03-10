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

function fullsys = couplefullsystem(Nb,Nt,Nsv, FILLING)
    params.Nb = Nb;
    params.Nt = Nt;
    params.Nsv = Nsv;
    params.FILLING = FILLING;
    fullsys.p = params;
    ebbeam = eulerbernoulli(Nb);
    tbeam = torsion(Nt);
    beam = phcat(ebbeam, tbeam);
    
    fullsys.ebbeam = ebbeam;
    fullsys.tbeam = tbeam;
    
    [pslosh p prigid] = dataexperiment(FILLING); % look dataexperiment.m file:
    fullsys.p.pslosh = pslosh;
    rbm = rb(prigid.mrb); %mass rigid body
    rbIb = rb(prigid.Irb+prigid.Ifluid);%bending rigid body
    rbIt = rb(prigid.Irb);%tosion rigid body
    rbc = phcat(rbIb,rbm);
    rbc = phcat(rbc,rbIt);
    fluid = rotsaintvenant(Nsv,FILLING);

    J = blkdiag(beam.J,fluid.J, rbc.J);
    Q =  blkdiag(beam.Q,fluid.Q, rbc.Q);
    JQ = J*Q;
    B = blkdiag(beam.B, fluid.B, rbc.B);
    D = blkdiag(beam.D, fluid.D, rbc.D);
    M = zeros(9,9); N = zeros(9,9);
    % u = [-mL      fL  |    -ML~T    | -FL, v0,  Mext  | MB      F     MT   ]
    % y = [-thetadL wdL |  -thedatL~T | -vL, F0, thetap | thetapB wB thetapT ]
    % u = [-e1^B(L)     -e1d^B(L)  |    e1(L)~T   |  e1(L)^F,  e2(0)^F,  Mext   | MB      F     MT   ]
    % y = [-e2d^B(L)     e2^B(L)   |  -e2(L)~T    | -e2(L)^F,  e1(0)^F, thetap  | thetapB wB thetapT ]
    M(1,8) = 1; M(1,2) = -1; % y8 = y2   
    M(2,8) = 1; N(2,5) = -1; % y8 = u5
    M(3,8) = 1; M(3,4) = 1; %  y8 = -y4

    M(4,7) = 1; M(4,1) = 1; % y7 = -y1  (bending angle RB and beam)

    M(5,9) = 1; M(5,3) = -1; % y9 = y3 (torsion angle RB and beam) 
    M(6,9) = 1; M(6,6) = -1; %  y9 = y6 (torsion angle RB and fluid)

    N(7,2) = 1; N(7,4) = -1; N(7,8) = 1; M(7,5) = 1; % sum of forces
    N(8,1) = -1; N(8,7) = 1; % sum of bending moments
    N(9,3) = 1; N(9,6) = 1; N(9,9) = 1; % sum of torsion moments


    A = [JQ, B; M*transpose(B)*Q, M*D+N];
    E = blkdiag(eye(size(JQ)), zeros(9,9));

    fullsys.A = A;
    fullsys.E = E;
    fullsys.M = M;
    fullsys.N = N;
    fullsys.Q = Q;
    fullsys.B = B;
    fullsys.D = D;
    fullsys.J = J;
end