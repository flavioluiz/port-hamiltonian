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

function exactfreqs = tb_exact(NUMTOR)
%% Calculate exact natural frequencies of torsion beam
% Function usage: exactfreqs = tb_exact(NUMTOR)
%   where NUMTOR is the number of torsion modes;
%   exactfreqs is a vector of frequencies in Hertz
%


% load experimental data
    [pslosh p prigid] = dataexperiment(0); % look dataexperiment.m file:
                                           % this loads all the
                                           % experimental data
                                           
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
    omegan = sqrt(GJ/Ip)/L*[1:2:(NUMTOR*2-1)]/2  /2;

     exactfreqs = omegan(:);
    
end