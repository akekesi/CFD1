%---------------------------------------
% CFD2_HA_03_Group L
%   XY            (xxxxxx)
%   XY            (xxxxxx)
%   Attila Kekesi (xxxxxx)
%
%   MATLAB R2020a
%---------------------------------------

% Lösung der Euler-Gleichungen auf nichtorthogonalen Gittern
% Simulation der Euler-Gleichungen in 2D, periodisch

%%
clear all
%close all
clc
clf
cla

vid = 0;                    % als Video zu speichern: vid = 1
vid_name = 'CFD2_HA03.avi'; % Name des Videos

%% Anfangswerte, Konstante
global params
alpha = 0.2;
beta  = 0.05;
delta = 0.2;
params.gamma = 1.4;
x0 = 2;
y0 = 2;
rho0 = 1.204;
p0 = 1e5;
u0 = 0;
v0 = 0;

N_xi = 100;
N_eta = 55;
xi_end = 2*pi;
eta_end = 2*pi;
h_xi = xi_end/N_xi;         % fuer Ableitungsmatrix
h_eta = eta_end/N_eta;      % fuer Ableitungsmatrix

xi = (1:N_xi)*2*pi/N_xi;
eta = (1:N_eta)*2*pi/N_eta;

x = @(xi,eta) xi-alpha*sin(xi-eta);
y = @(xi,eta) eta-alpha*sin(xi-eta);
rho = @(x,y) rho0 + beta*exp(-((x-x0).^2 + (y-y0).^2)/delta);
p = @(rho) p0*(rho/rho0).^params.gamma; % Isentropen-Beziehung
c = @(p,ro) sqrt(params.gamma*p./ro);   % Schallgeschwindigkeit

dt = 0.0002;    % Zeitschritt
N = 150;        % Anzahl der Zeitschritte

%% (a1) Erzeugung eines periodischen Gitters in xi und eta
[XI,ETA] = meshgrid(xi,eta);
X = x(XI,ETA);
Y = y(XI,ETA);

%% (a2) Berechnung der Geometriefaktoren
D_TW_xi = D_TW(N_xi,h_xi);
D_TW_eta = D_TW(N_eta,h_eta);
params.D_kron_xi = kron(speye(N_eta),D_TW_xi);
params.D_kron_eta = kron(D_TW_eta,speye(N_xi));

params.ex_xi = 1+D_xi(X-XI);
params.ey_xi = D_xi(Y);
params.ex_eta = D_eta(X);
params.ey_eta = 1+D_eta(Y-ETA);

params.g11 = params.ex_xi.^2 + params.ey_xi.^2;                         % 11:       xixi
params.g22 = params.ex_eta.^2 + params.ey_eta.^2;                       % 22:       etaeta
params.g12 = params.ex_xi.*params.ex_eta + params.ey_xi.*params.ey_eta; % 12=21:    xieta = xieta
params.J = params.ex_xi.*params.ey_eta - params.ex_eta.*params.ey_xi;

%% (c) Berechnung der Anfangsbedingung
Rho0 = rho(X,Y);
P0 = p(Rho0);
C0 = c(P0,Rho0);
U0 = u0*ones(size(X));
V0 = v0*ones(size(X));
q0(:,:,1) = sqrt(Rho0);
q0(:,:,2) = sqrt(Rho0).*U0;
q0(:,:,3) = sqrt(Rho0).*V0;
q0(:,:,4) = P0;

%% (d) Bestimmung der CFL-Zahl
% ???
%CFLx = u*dt/dx;
%CFLy = v*dt./dy;
%CFLc = c*dt./min(h_xi,h_eta)
% ???

%% Plot
figure(1)
pcolor(X,Y,C0)
colorbar
colormap jet
title('Schallgeschwingikeit - Anfangswerte')
simulation(X,Y,q0,dt,N,vid,vid_name);

%% (b) Rechte Seite
function rhs = rhs_ns2d(q)
global params
rho = q(:,:,1).^2;
u = q(:,:,2)./q(:,:,1);
v = q(:,:,3)./q(:,:,1);
p = q(:,:,4);
u_ = u.*params.ey_eta - v.*params.ex_eta;
v_ = v.*params.ex_xi - u.*params.ey_xi;

rhs = zeros(size(q));
rhs(:,:,1) = -(D_xi(u_.*rho) + D_eta(v_.*rho))./(2*params.J.*sqrt(rho));
rhs(:,:,2) = -(1/2*(D_xi(u_.*rho.*u) + rho.*u_.*D_xi(u) + D_eta(v_.*rho.*u) + rho.*v_.*D_eta(u))+...
    D_xi(params.ey_eta.*p) - D_eta(params.ey_xi.*p))./(params.J.*sqrt(rho));
rhs(:,:,3) = -(1/2*(D_xi(u_.*rho.*v) + rho.*u_.*D_xi(v) + D_eta(v_.*rho.*v) + rho.*v_.*D_eta(v))+...
    D_eta(params.ex_xi.*p) - D_xi(params.ex_eta.*p))./(params.J.*sqrt(rho));
rhs(:,:,4) = -(params.gamma*(D_xi(u_.*p) + D_eta(v_.*p)) -...
    (params.gamma-1)*u.*D_xi(params.ey_eta.*p) - u.*D_eta(params.ey_xi.*p) + v.*D_eta(params.ex_xi.*p) - v.*D_xi(params.ex_eta.*p))./params.J;
end

%% (e) Simulation
function [] = simulation(X,Y,q0,dt,N,vid,vid_name)
% Video
if vid == 1
    Video = VideoWriter(vid_name);
    Video.FrameRate = 13;
    open(Video)
end
for i = 1:1:N+1
    if i == 1
        q = q0;
    else
        q = RK4(q,dt,@rhs_ns2d);
    end
    rho = q(:,:,1).^2;
    u = q(:,:,2)./q(:,:,1);
    v = q(:,:,3)./q(:,:,1);
    p = q(:,:,4);

    figure(2)
    sgtitle(['t = ',num2str((i-1)*dt),' s'])

    subplot(2,2,1)
    pcolor(X,Y,rho)
    title("Dichte [kg/m^3]")
    axis square
    colorbar

    subplot(2,2,2)
    pcolor(X,Y,p)
    title("p [Pa]")
    axis square
    colorbar

    subplot(2,2,3)
    pcolor(X,Y,u)
    title("u [m/s]")
    axis square
    colorbar

    subplot(2,2,4)
    pcolor(X,Y,v)
    title("v [m/s]")
    axis square
    colorbar

    drawnow
    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end
if vid == 1
    close(Video)
end
end

%% Ableitungen fuer beliebige Funktionen/Felder F
function F_xi = D_xi(F)
global params
D_kron_xi = params.D_kron_xi;
[n,m] = size(F);
F = reshape(F.',n*m,1);
F_xi = D_kron_xi*F;
F_xi = reshape(F_xi,m,n)';
end

function F_eta = D_eta(F)
global params
D_kron_eta = params.D_kron_eta;
[n,m] = size(F);
F = reshape(F.',n*m,1);
F_eta = D_kron_eta*F;
F_eta = reshape(F_eta,m,n)';
end

%% Expicit Derivate Matrix with improved wave numbers of Tam and Webb
% from isis
function D = D_TW(N,h)
% Expicit Derivate Matrix with improved wave numbers of Tam and Webb
% Periodic Version
% D = D_TW(N,h)
% N number of Points
% h = Delta x   

alpha2=-1/(h)*0.18941; %*(2*pi)/h;
alpha3= 1/(h)*0.02652; % *(2*pi)/h;
alpha1= 1/(2*h) - 2*alpha2 - 3*alpha3;

D= alpha1*(diag(ones(N-1,1),1)-  diag(ones(N-1,1),-1)  + diag(ones(1,1),-N+1) -  diag(ones(1,1),N-1))    +...  
   alpha2*(diag(ones(N-2,1),2)-  diag(ones(N-2,1),-2) +  diag(ones(2,1),-(N-2))- diag(ones(2,1),N-2))  +...        
   alpha3*(diag(ones(N-3,1),3)-  diag(ones(N-3,1),-3) +  diag(ones(3,1),-(N-3))- diag(ones(3,1),N-3) );

%D= 1/(h)*D ; 
end

%% Runge Kutta 4. Ordnung
% from isis
function u1 = RK4(u0,dt,rhs)
% Runge Kutta 4. Ordnung
% die vier Aufrufe der rechten Seiten; 
% egal, ob rhs zahl, vektor  oder vektoren=matrix zurueckgibt: 
% solange Dimension u entspricht geht alles glatt.  
 
k1=rhs(u0        );
k2=rhs(u0+dt/2*k1);
k3=rhs(u0+dt/2*k2);
k4=rhs(u0+dt  *k3);

u1 = u0+dt/6*(k1+2*k2+2*k3+k4);
end