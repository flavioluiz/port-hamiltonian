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

%% this file show some results of the full coupled system
%
% This file is separated in three main blocks:
%   1) load the coupled system and calculate
%       eigenvalues;
%   2) make animations of the eigenvectors;
%   3) simulate the system and make animation;
%
%  once part 1 is loaded, user can try part 2 and 3
%
%  I should improve the usage of this code later...

%% PART 1  -- LOADING AND EIGENVALUES CALCULATION
% This first block is used to load the coupled system from 
% function coupledfullsystem and finding the system eigenvalues
% and eigenvectors.
clear all;
addpath('.\main');


Nb = 20; % number of bending modes
Nt = 20; % number of torsion modes
Nsv = 20; % number of sloshing modes
%ALERT: current version of 3D plot function needs
%       to use Nb = Nt (same number of elements
%                       for torsion and bending)
FILLING = 0.75; % tank filling ratio

% load the system
fullsys = couplefullsystem(Nb,Nt,Nsv, FILLING);

A = fullsys.A;
E = fullsys.E;

% take eigenvectors and eigenvalues
[autovec autoval] = eig(A,E);
autoval = diag(autoval);
% only non-infinity eigenvalues are taken
%     *inf valued eigenvalues are related with constraints
noninfindex = abs(real(autoval))<inf;
noninfautoval = autoval(noninfindex);
noninfautovec = autovec(:,noninfindex);
[autovalsort ord] = sort(abs(imag(noninfautoval)));
autovecsort = noninfautovec(:,ord);
% show first natural frequencies, in Hertz
(autovalsort(6:2:24))'/2/pi
    % notice that eigenvalues from 1 to 5
    % are zero: they are related with rigid body
    % integration (states that integrate a speed
    % to find a position)

% check max real part of eigenvalues:
max(real(autoval(real(autoval)<inf)))
   % if everything is ok, this should be ZERO

%% PART 2
%% this saves a video of mode ii in mp4
ii = 2;  % number of mode
freqn = autovalsort(ii)/2/pi
SV1 = (2*Nb+2*Nt+1):(2*Nb+2*Nt+Nsv);
writerObj = VideoWriter('results/mode2.mp4','MPEG-4');
writerObj.Quality = 100;
open(writerObj);
figure(1);
autovecn = autovecsort(:,ii)/max(abs(autovecsort(SV1,ii)));
tt = 0;
plotfull(0.004*(real(autovecn)*cos(2*pi*tt) +  imag(autovecn)*sin(2*pi*tt) ),fullsys.p); 
axis tight
opengl('software')
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','opengl');
set(gcf, 'OuterPosition',  [846         122        1073         959]);
for tt = 0:(1/30):(2/freqn)
    figure(1);
    plotfull(0.004*(real(autovecn)*cos(2*pi*freqn*tt) +  imag(autovecn)*sin(2*pi*freqn*tt) ),fullsys.p);    
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);

%% this saves a video of mode ii in GIF
ii = 12; % number of mode
freqn = autovalsort(ii)/2/pi
SV1 = (2*Nb+2*Nt+1):(2*Nb+2*Nt+Nsv);
filename = 'results\mode3.gif';

figure(1);
autovecn = autovecsort(:,ii)/max(abs(autovecsort(SV1,ii)));
tt = 0;
set(gcf, 'OuterPosition',  [1211         190         670         360]);
for tt = 0:(1/30):(4/freqn)
    figure(1);
    plotfull(0.004*(real(autovecn)*cos(2*pi*freqn*tt) +  imag(autovecn)*sin(2*pi*freqn*tt) ),fullsys.p);    
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,128);
    if tt == 0;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',1/30);
    end
end

%% PART 3
%% Simulation:
idyn =@(t,x,xp) (E*xp-A*x);

x0 = zeros(size(A,1),1);
% defining initial condition:
%  stationary case, with 4 Newtons tip force is considered:
x0b = - (fullsys.ebbeam.J*fullsys.ebbeam.Q) \ fullsys.ebbeam.B(:,2) * 4;
x0(1:Nb*2) = x0b;
[T, X] = ode15i(idyn, [0:0.05:10], x0,x0*0);

%% video/ animated GIF of simulation results
HH = [];
filename = 'results\simulation.gif';
figure(1);
set(gcf, 'OuterPosition',  [1211         190         670         360]);

for i = 1:1:length(T)    
    plotfull(X(i,:)',fullsys.p);

    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',0.05);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',0.05);
    end
end
%% calculate system energy variation
for i = 1:length(T)
    HH(i) = X(i,1:length(fullsys.Q))*(fullsys.Q)*X(i,1:length(fullsys.Q))';
end
figure; plot(T,(HH-HH(1))/HH(1)*100)