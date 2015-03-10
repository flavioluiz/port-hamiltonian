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

function exactfreqs = exact_slosh(NUM,FILLING)
%% Calculate exact natural frequencies of shallow water equations
% Function usage: exactfreqs = exact_slosh(NUM,FILLING)
%   where NUM is the number of sloshing modes;
%       FILLING is the tank filling ratio (0 to 1)
%   exactfreqs is a vector of frequencies in Hertz
%


% load experimental data
    [pslosh p prigid] = dataexperiment(FILLING); % look dataexperiment.m file:
                                           % this loads all the
                                           % experimental data
    
    omegan = pi*sqrt(pslosh.g*pslosh.h)/pslosh.a*[1:NUM]/2/pi;

     exactfreqs = omegan(:);
    
end