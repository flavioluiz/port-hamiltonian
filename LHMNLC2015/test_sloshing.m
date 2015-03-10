%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%																   %
%   Flavio Luiz Cardoso-Ribeiro: http://github.com/flavioluiz/     %
%			ISAE-Supaero   / Instituto Tecnologico de Aeronautica  %
%								CNPq   - Brazil 				   %
%																   %
%    This project is part of ANR Project HAMECMOPSYS:              %
%					https://hamecmopsys.ens2m.fr/ 				   %
%																   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sloshing test file:
%  The following script tests rotsaintvenant.m file, by comparing
%  numerically obtained eigenvalues against known exact solution
%
clear all;
FILLING = 0.4;
addpath('.\main');
i = 1;
NUMELEMVEC = [10 20 40 100];
for NUMELE = NUMELEMVEC
    % load beam
    fluid = rotsaintvenant(NUMELE,FILLING);
	M = [1, 0, 0; 0, 0, 0; 0, 0, 1];
	N = [0, 0, 0; 0, 1, 0; 0, 0, 0];
    A = [fluid.J*fluid.Q, fluid.B;
			M*fluid.B'*fluid.Q,M*fluid.D+N];
	E = eye(NUMELE*2+2+3);
	E(end-2:end,end-2:end) = 0;
	
    natfreqs = sort(damp(eig(A,E))/2/pi);     
    % 
    freqs(:,i) = natfreqs(3:2:16);
    i = i+1;
end

% print results
values = [freqs exact_slosh(7,FILLING)];
fprintf(' Number of elements: | %d | %d | %d | %d | Exact \n', NUMELEMVEC(1), NUMELEMVEC(2),...
                                    NUMELEMVEC(3), NUMELEMVEC(4));
fprintf(' --- | --- | --- | --- | --- | --- \n');
for iline = 1:7
    fprintf('   | %.4f | %.4f | %.4f | %.4f | %.4f \n', values(iline,1) , values(iline,2), ...
                        values(iline,3), values(iline,4), values(iline,5));   
end