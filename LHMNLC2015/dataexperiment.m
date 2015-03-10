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

function [pslosh pbeam prigid] = dataexperiment(filling)
    %% internal tank ("sloshing related" parameters):
    % [pslosh pbeam prigid] = dataexperiment(filling)
    pslosh.g = 9.8; % m/s^2
    pslosh.rho = 999; % water density: kg/m^3      %    +- 1
    L = 0.470;    % tank internal length (m)       %   +- 0.001
    R = 0.105/2; % tank internal radius (m)        %   +- 0.0002
    if nargin < 1
    e = 0.7;     % tank filling ratio (in height: h/2R)
    else
        e = filling;
    end
    %"equivalent" rectangular tank (same surface area and volume):
    theta = 2*acos(1-2*e);
    ls = 2*R*sqrt(1-(2*e-1)^2);
    Vcyl = R^2/2*(theta-sin(theta))*L;

    pslosh.a = L;           % rectangular tank length
    pslosh.b = ls;          % rectangular tank width
    pslosh.h = Vcyl/ls/L;   % rectangular tank liquid height
    
    pslosh.R = R;
    pslosh.H = e*2*R;
    %% structural parameters
    %plate parameters   - data in Ref. Robu2012,  'Simultaneous Hinf Vibration
    % Control of Fluid/Plate via Reduced-Order Controller')
    % --- IEEE Transactions on Control Sys Tech:
    E = 75e9; %Young modulus
    mu = 0.33;  %Poisson ratio
    rho = 2970; %density
    L = 1.36; %length
    h = 0.005;  %thickness       +- 0.00005
    l = 0.16;   %width           +- 0.001

    mrb = 0; % rigid body (tank) mass (here it is zero since we are adding
    % all rigid bodies when coupling with sloshing ("D" matrix))
    Irb = 0; % rigid body (tank) inertia (along axis x and y)
    xrb = 1.28; % position of rigid body (along axis x)        +- 0.001

    pbeam.E = E; pbeam.mu = mu; pbeam.L = L; pbeam.rho = rho; pbeam.h = h;
    pbeam.l = l; pbeam.mrb = mrb; pbeam.Irb = Irb; pbeam.xrb = xrb;

    % piezoelectric data (data in ref. Robu 2012)
    rhop = 7800; % density
    Ep = 67e9; % young modulus
    tp = 5e-4; % thickness
    d31 = -2.1e-10; % piezoelectric constant
    pbeam.Ep = Ep; pbeam.tp = tp; pbeam.d31 = d31; pbeam.rhop = rhop;

    %% rigid mass and inertias
    %(tank, fixation rings and fluid inertia in "y" direction)    
    %tank:
    
    % tank data - data in Robu2012/IEEE
    rhoplast = 1180; % plastic density (kg/m^3)
    Der = 0.11;      % tank external diameter (m)
    Dir = 0.105;     % tank internal diameter (m)
    Lr = 0.50;       % tank length (m)              +- 0.001

    % fixation rings data - data measured/approximated
    % obs.: external diameter and thickness doesn't represent the exact
    % values. actually, the fixation rings have several holes/screws.
    % Def, Di and hf are approximate values. They were adjusted so that
    % the TOTAL fixation ring masses are the measured ones
    %   (588.6g)
    rhofix = 2970;       % fixation rings density (= plate,kg/m^3)
    Def = 0.1445 + 0.15*0; % fixation rings external diameter (m)
    Dif = 0.11;            % fixation rings internal diameter (m)
    hf = 0.02876+0.0325*0; % fixation rings thickness (m)
    distf = (h+hf)/2;      % fixation rings position (relative to plate mean axis) (m)

    Dh = Der;        % PLATE hole diameter (m)
    massrlat = rhoplast* pi * (Der^2-Dir^2)/4 * Lr; % tank lateral mass 
    massrtip = 0.313/2;     % each tip mass (measured)
                            % OBS.: the tip mass was obtained by measuring
                            % the tank mass and comparing with the predicted
                            % lateral tank mass. In order to confirm if the
                            % lateral tank mass was correct, I also
                            % compared the results with a 30cm tank (in
                            % addition to the 50cm one).
                            % The tip consists of: plastic disc + glue +
                            % plug: not so simple geometry, but the effects
                            % on the system mass/inertia can be
                            % approximated as a point mass.
    massfix = rhofix*pi*(Def^2-Dif^2)/4 * hf;         % each fixation rings total
    masshole = rho*pi*Dh^2/4 * h;                     % plate hole mass - it will reduce the rigid unit mass!
    mass = (massrlat + massrtip*2+2*massfix-masshole);% total rigid added mass (kg)
    
    % calculate rigid body inertia (x and y direction, since tank and
    % fixation rings are symmetric):
    Inertialat = (massrlat/12*(3*(Der^2+Dir^2)/4+Lr^2)); % tank lateral inertia
    Inertiatip = massrtip*(Lr/2-0.01)^2; % tank tip inertia
    Inertiafix = massfix*(3*(Dif^2+Def^2)/4+hf^2)/12+massfix*distf^2; % fixation rings inertia
    Inertiahole = masshole*(Dh^2)/16;   % inertia of hole in the plate
    Inertia = (Inertialat + 2*Inertiatip + 2*Inertiafix - Inertiahole); %(kg.m^2)
 
    % "frozen" fluid inertia (y direction)
    massfluid = pslosh.a*pslosh.b*pslosh.h*pslosh.rho;
    IFLUID = massfluid*pslosh.a^2/12;    %(kg.m^2)
    
    prigid.xrb = xrb; %position of rigid units (m)
    prigid.mrb = mass; %tank + fixation rings mass (kg)
    prigid.Irb = Inertia; %tank + fixation rings inertia ("y" and "z" directions) (kg.m^2)
    prigid.Ifluid = IFLUID; % "frozen" fluid inertia in "y" direction (kg.m^2)
    prigid.mfluid = massfluid; % "frozen" fluid mass
    
    %% auxiliary equations and parameters: most equations are implemented on
    % other codes
    % they are repeated here only for verification purposes

    % beam structural parameters
    J = 0.33*l*h^3;              % torsion constant
    IY = rho/12*(l^2 + h^2)*l*h; % section moment of inertia (X direction)  % VERIFICAR!!!
    Iz = l*h^3/12;               % section area inertia (X direction)
    Ix = l^3*h/12;               % section area inertia (X direction)
    G = E/2/(1+mu);              % shear modulus
   
end

%% measures:
% large tank: 50 cm  (47 cm liquid)
% tank filled up to 34.3 cm (+-0.2):
%    mass = 3820.3 g +- 0.1 g
% tank filled 100 %:
%    mass = 4938.0 g +- 0.4 g
% tank filled 0%:
%    mass = 821.8 g +- 0.1 g
%
% small tank: 30 cm (27 cm liquid)
%  empty tank:
%    606.7 g +- 0.1 g
%
% small tank: 30 cm (27 cm liquid): broken one
%    578.1 g +- 0.1 g
%

%% new measures:
% TANK WEIGHT:
%   * empty: 818.2 gramas (erro do aparelho: 0.1, mas gotas extras
%   podem adicionar uns 0.2 g)
%   * 0 +-0.1cm:      839.8  +-0.1g
%   * 20.2 +- 0.2cm: 2551.2  +-0.2g
%   * 35.0 +- 0.1cm: 3879.2  +-0.2g
%   * 47.0 +- 0.1cm: 4933.4  +-0.2g
%   * full: 4939.9 +- 0.1 g (mas nao estava 100% cheio, impossivel
%   fazer isso)
%
% TANK medidas dos bouchons individuais: 15.9, 15.8, 16.7, 16.2
% (gramas)
%
% TANK dimensoes:
% L = 30.0 +-0.1 cm
% Li = 27.0 +- 0.1cm
% espessura do disco interno de plastico: 0.70 +-0.05 cm
% distancia entre a parte interna do disco e o fim do bouchon: 4.5cm
% distancia entre a parte interna do disco e a extremidade: 1.5cm
%
%
% ANEIS de fixacao:
%  1) lado direito (visto de frente)
%       561.4 g (sem parafusos)
%       568.2 g (parafuso central)
%       580.6 g (parafuso central+interm)
%       596.2 g (parafuso central+interm+sup)
%   2) lado esquerdo
%       565.3 g (sem parafusos)
%       580.9 g (parafuso superior)
%
%
% peso dos parafusos do anel:
%   central: 3.3g
%   interm: 6.2g
%   sup: 15.7g
%
% dimensoes do anel:
%  diametro exterior 1: 16.0cm, 16.0cm, 16.0cm +- 0.1
%  diametro exterior 2 (anel 1): 145.46mm, 145.52mm,145.45mm,145.40mm
%  diametro exterior 2 (anel 2): 144.62mm, 144.56mm,144.51mm
%  diametro interior (anel 1): 110.53mm, 110.8mm,110.81mm,110.72mm
%  diametro interior (anel 2): 110.73mm, 110.90mm,110.62mm
%  altura placa fina: 5.16 5.22 5.02 4.97 4.48mm
%  altura total anel: 30.15 30.01 30.13 29.79mm
%
%
% DIMENSOES DA PLACA
%    comprimento: 136.3 +-0.2 cm
%    largura: 16.0 +- 0.1 cm
%    espessura: 4.95mm, 4.98mm, 4.92mm, 4.92mm, 4.84mm, 5.02mm, 4.88mm,
%    4.86mm
%
%
% POSICAO DO FURO NA PLACA:
%  relativo a ponta da placa:
%    entre 2.4cm e 13.65cm
%  diametro interno do furo:
%          112.05 +-0.3mm



