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

function sys = phcat(sys1,sys2)
 % concatenate port hamiltonian systems
 sys.J = blkdiag(sys1.J,sys2.J);
 sys.Q = blkdiag(sys1.Q,sys2.Q);
 sys.B = blkdiag(sys1.B,sys2.B);
 sys.D = blkdiag(sys1.D,sys2.D);
end