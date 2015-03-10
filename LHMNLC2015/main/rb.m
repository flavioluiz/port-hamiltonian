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

function rb = rb(mass)
    % rigid body writen in port Hamiltonian framework
    % input -> force(moment)
    % output-> speed(angle rate)
    % Xp = JQX + Bu
    % y = B'QX
    rb.J = [0,-1;1,0];
    rb.Q = diag([1/mass,0]);
    rb.B = [1;0];
    rb.D = 0;
end