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

function plotfull(X,params)
%% PLOT full coupled system
% This file plot the full coupled system, given 
%                                     X:  the state vector;
%                                     params: params structure, from
%                                             couplefullsystem object
%
%
    Nb = params.Nb;
    Nt = params.Nt;
    Nsv = params.Nsv;
    FILLING = params.FILLING;
    pslosh = params.pslosh;
    EB1 = 1:(Nb);
    EB2 = (Nb+1):(2*Nb);
    TB1 = (2*Nb+1):(2*Nb+Nt);
    TB2 = (2*Nb+Nt+1):(2*Nb+2*Nt);
    SV1 = (2*Nb+2*Nt+1):(2*Nb+2*Nt+Nsv);
    SV2 = (2*Nb+2*Nt+Nsv+1):(2*Nb+2*Nt+2*Nsv);
    RB = (2*Nb+2*Nt+2*Nsv+1):(2*Nb+2*Nt+2*Nsv+6);
    DZB = 1.36/Nb;
    DZF = 0.47/Nsv;
    
    
%     subplot(2,1,1);
    W = finddisplacements(X(EB1)/DZB);
    TO = findtorsionangle(X(TB1)/DZB);
    x = linspace(0,1.36,length(TO));
    y = linspace(-0.08,0.08,10);
    WW = (TO*y + repmat(W,1,10));
    [XX YY] = meshgrid(x,y);
    tc = 0:pi/10:2*pi;
    
    [Xc,Yc,Zc] = cylinder(cos(tc)*0+pslosh.R);
    Zc = Zc*0.47-0.47/2;
    
    Xc = Xc + 1.36;
    Zc = Zc + W(end);
    mesh(Xc,Zc, Yc,'FaceColor','red','EdgeColor','none'); alpha(0.2); hold on;
    %camlight left;
    %lighting phong

    mesh(XX,WW',YY);    
    plotar = SV1;
    
    [XXsl YYsl] = meshgrid(linspace(-0.47/2,0.47/2,length(plotar)),linspace(-pslosh.b/2,pslosh.b/2,4));
    WWsl = repmat(X(SV1)'/DZF*pslosh.h,4,1)+pslosh.H-pslosh.R;
    YYsl = YYsl+1.36;
    XXsl = XXsl+W(end);
    surf(YYsl,XXsl,WWsl,'FaceColor','blue','EdgeColor','none');
    %oldpos = get(gca,'CameraPosition')
    

    newpos = [8.6157   -3.0749    1.5935];%positioned at z=10
    target = [0.5 0 0];%aimed at the origin
    set(gca,'CameraPosition',newpos,'CameraTarget',target);
    hold off;
    axis equal
    axis([0 1.36 -0.53 0.53 -0.08 0.08])

end

function w = finddisplacements(d2wdx2)
% this function finds the displacement vector w(z,t), from d2wdz2
    N = length(d2wdx2);
    dz = 1.36/N;
    w = dz*tril(ones(N),-1)*integrated2wdz2(d2wdx2) + dz^2/2 * tril(ones(N))*d2wdx2;
end
function dwdz = integrated2wdz2(X)
% this function calculates dwdz from d2wdz2
% it basically integrate d2wdz2, which is constant in each element
    dwdz = 1.36/length(X)*tril(ones(length(X)))*X;
end

function theta = findtorsionangle(dthetadz)
% this function finds the torsion angle from dthetadz vector
% it integrates dthetadz, which is constant in each element
    theta = 1.36/length(dthetadz)*tril(ones(length(dthetadz)))*dthetadz;
end