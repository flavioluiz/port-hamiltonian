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

%% Euler-Bernoulli beam coupled with rigid body test file:
%  The following script tests eulerbernoulli.m file, by comparing
%  numerically obtained eigenvalues against ______
%
clear all;
addpath('.\main');
i = 1;

M = zeros(4,4); N = zeros(4,4);
N(1,1) = -1; N(1,3) = 1; % sum of moments ueb(1), urbB(1)
M(2,1) = -1; M(2,3) = -1; % thetadL = thedadRB
M(3,2) = 1; M(3,4) = -1; % wL = Wrb
N(4,2) = 1; N(4,4) = 1; % sum of forces

NUMELEMVEC = [10 20 40 100];
for NUMELE = NUMELEMVEC
    % load beam
    eb = eulerbernoulli(NUMELE);
    rbIb = rb(0.0319);%rotation inertia rigid body
    rbm = rb(1.8481); %mass rigid body
    rbc = phcat(rbIb, rbm);
    
    sys = phcat(eb,rbc);
    
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
values = [freqs exact_eb(7)];
fprintf(' Number of elements: | %d | %d | %d | %d \n', NUMELEMVEC(1), NUMELEMVEC(2),...
                                    NUMELEMVEC(3), NUMELEMVEC(4));
fprintf(' --- | --- | --- | --- | --- | --- \n');
for iline = 1:7
    fprintf('   | %.4f | %.4f | %.4f | %.4f \n', values(iline,1) , values(iline,2), ...
                        values(iline,3), values(iline,4));   
end