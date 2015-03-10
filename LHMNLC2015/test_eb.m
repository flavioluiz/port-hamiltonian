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

%% Euler-Bernoulli beam test file:
%  The following script tests eulerbernoulli.m file, by comparing
%  numerically obtained eigenvalues against known exact solution
%
clear all;
addpath('.\main');
i = 1;
NUMELEMVEC = [10 20 40 100];
for NUMELE = NUMELEMVEC
    % load beam
    eb = eulerbernoulli(NUMELE);
    % find eigenvalues of fixed-free beam:
    natfreqs = sort(damp(eb.J*eb.Q)/2/pi);     
    % 
    freqs(:,i) = natfreqs(1:2:14);
    i = i+1;
end

% print results
values = [freqs exact_eb(7)];
fprintf(' Number of elements: | %d | %d | %d | %d | Exact \n', NUMELEMVEC(1), NUMELEMVEC(2),...
                                    NUMELEMVEC(3), NUMELEMVEC(4));
fprintf(' --- | --- | --- | --- | --- | --- \n');
for iline = 1:7
    fprintf('   | %.4f | %.4f | %.4f | %.4f | %.4f \n', values(iline,1) , values(iline,2), ...
                        values(iline,3), values(iline,4), values(iline,5));   
end