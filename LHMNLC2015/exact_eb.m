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

function exactfreqs = eb_exact(NUMBEN)
%% Calculate exact natural frequencies of euler-bernoulli beam
% Function usage: exactfreqs = eb_exact(NUNBEN)
%   where NUNBEN is the number of bending modes;
%   exactfreqs is a vector of frequencies in Hertz
%


% load experimental data
    [pslosh p prigid] = dataexperiment(0); % look dataexperiment.m file:
                                           % this loads all the
                                           % experimental data
                                           
    E = p.E; % Aluminium Young Modulus (Elasticity modulus): N/m^2
    mu = p.mu; % Poisson's ratio (for the aluminium)
    L = p.L; % beam Length: m
    rho = p.rho;  % aluminium density: kg/m^3
    h = p.h; % beam thickness: m
    l = p.l; % beam width: m

    S = h*l;  % beam sectional area m^2
    I = l*h^3/12; % sectional area inertia

    Ep = p.Ep; % piezoelectric Young modulus
    tp = p.tp; % piezoelectric patch thickness
    d31 = p.d31; % piezoelectric constant

    mup = mu; % piezoelectric patch
    
    
    %modal frequencies omegan
        betan0 = [0.1 3 6 8 10 13 15 17 19.5 22 24 26 29 31 33.3 35.6];
                % betan0: starting values for the fsolve algorithm
        for i = 1:NUMBEN
            betan(i) = fsolve(@(x) cos(x*L)*cosh(x*L) + 1,betan0(i));            
        end
        omegan = (betan).^2 * sqrt(E*I/(rho*S)); % this finds a VECTOR of omegan (natural frequencies)
    
     exactfreqs = omegan'/2/pi;
    
end