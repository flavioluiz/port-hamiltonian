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

%% Torsion beam coupled with rigid body test file:
%  The following script tests eulerbernoulli.m file, by comparing
%  numerically obtained eigenvalues against ______
%
clear all;
addpath('.\main');
i = 1;

M = zeros(2,2); N = zeros(2,2);
M(1,1) = 1; M(1,2) = 1; % thetadRB, -thetadT_L
N(2,1) = 1; N(2,2) = -1; % sum of torque: m^RB, m^T_L

NUMELEMVEC = [10 20 40 100];
for NUMELE = NUMELEMVEC
    % load beam
    tb = torsion(NUMELE);
    rbIt = rb(0.0319);%tosion rigid body
    
    sys = phcat(tb,rbIt);
    
    JQ = sys.J*sys.Q;
    Q = sys.Q;
    B = sys.B;
    D = sys.D;
    A = [JQ, B; M*transpose(B)*Q, M*D+N];
    E = blkdiag(eye(size(JQ)), zeros(size(M)));
    %sort(damp(eig(A,E))/2/pi)


    %
    [autovec autoval] = eig(A,E);
    autoval = diag(autoval);
    %real(autoval)
    noninfindex = abs(real(autoval))<inf;
    noninfautoval = autoval(noninfindex);
    noninfautovec = autovec(:,noninfindex);
    [autovalsort ord] = sort(abs(imag(noninfautoval)));
    autovecsort = noninfautovec(:,ord);
    freqs(:,i) = (autovalsort(3:2:16))'/2/pi

    i = i+1;
end

% print results
values = [freqs];
fprintf(' Number of elements: | %d | %d | %d | %d \n', NUMELEMVEC(1), NUMELEMVEC(2),...
                                    NUMELEMVEC(3), NUMELEMVEC(4));
fprintf(' --- | --- | --- | --- | --- | --- \n');
for iline = 1:7
    fprintf('   | %.4f | %.4f | %.4f | %.4f \n', values(iline,1) , values(iline,2), ...
                        values(iline,3), values(iline,4));   
end