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

function sys = massspring(mass,K)
    % mass spring written in port Hamiltonian framework
    % input -> force(moment)
    % output-> speed(angle rate)
    % Xp = JQX + Bu
    % y = B'QX
    sys.J = [0,-1;1,0];
    sys.Q = diag([1/mass,K]);
    sys.B = [1;0];
    sys.D = 0;
end