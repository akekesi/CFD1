%---------------------------------------
% CFD2_HA_02_Group L
%   XY            (xxxxxx)
%   XY            (xxxxxx)
%   Attila Kekesi (xxxxxx)
%
%   MATLAB R2020a
%---------------------------------------

% Lösung einfacher Gleichungen auf nichtorthogonalen Gittern
% Simulation eines passiven Skalars in 2D

%%
clear all
%close all
clc
clf
cla

%% Anfangswerte, Konstante
global params
alpha = 0.2;
N_xi = 100;                     % ??? Wenn es klein ist --> nicht schoen :) ???
N_eta = 55;                     % ??? Wenn es klein ist --> nicht schoen :) ???
xi_end = 2*pi;
eta_end = 2*pi;
h_xi = xi_end/(N_xi);         % fuer Ableitungsmatrix
h_eta = eta_end/(N_eta);      % fuer Ableitungsmatrix

xi = (1:N_xi)*2*pi/N_xi;
eta = (1:N_eta)*2*pi/N_eta;

x = @(xi,eta) xi-alpha*sin(xi-eta);
y = @(xi,eta) eta-alpha*sin(xi-eta);

params.u = 1;
params.v = 1;
params.gamma = 0.002;
dt = 0.05;   % Zeitschritt
N = 150;     % Anzahl der Zeitschritte

%% (a) Erzeugung eines periodischen Gitters in xi und eta
[XI,ETA] = meshgrid(xi,eta);
X = x(XI,ETA);
Y = y(XI,ETA);

%params.u = u.*ones(size(X));    % ??? es ist jetzt nicht noetig, da u ueberall gleich ist ???
%params.v = v.*ones(size(X));    % ??? es ist jetzt nicht noetig, da v ueberall gleich ist ???
%% (b) Berechnung der Geometriefaktoren
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

%% Plot
[X,Y,PHI] = simulation(X,Y,dt,N);
[~,~,N] = size(PHI);
figure(1)
for i = 1:1:N
    pcolor(X,Y,PHI(:,:,i))
    xlabel('x')
    yl = ylabel('y');
    set(yl,'rotation',0)
    title('Simulation eines passiven Skalars in 2D')
    colormap jet
    drawnow
end

%% (c) Transportoperator
function uGradPhi = u_grad_phi(u,v,phi)
global params
phi_xi = D_xi(phi);
phi_eta = D_eta(phi);

gradPhi1 = params.ey_eta.*phi_xi;
gradPhi2 = -params.ex_eta.*phi_xi;
gradPhi3 = -params.ey_xi.*phi_eta;
gradPhi4 = params.ex_xi.*phi_eta;

u_eff = u.*params.ey_eta-v.*params.ex_eta;
v_eff = -u.*params.ey_xi+v.*params.ex_xi;

uGradPhi = 1./params.J.*(u_eff.*(gradPhi1+gradPhi3)+v_eff.*(gradPhi2+gradPhi4));
end

%% (d) Laplace Operator
function LaplacePhi = Laplace(phi)
global params
phi_xi = D_xi(phi);
phi_eta = D_eta(phi);
LaplacePhi = 1./params.J.*...
    (D_xi(1./params.J.*(params.g22.*phi_xi-params.g12.*phi_eta))+...
    D_eta(1./params.J.*(-params.g12.*phi_xi+params.g11.*phi_eta)));
end

%% (e) Rechte Seite der Transportgleichung
function rhsPhi = rhs_PasScalar(phi)
global params
rhsPhi = -u_grad_phi(params.u,params.v,phi)+params.gamma*Laplace(phi);
end

%% (f) Simulation der Gleichung
function [X,Y,PHI] = simulation(X,Y,dt,N)
x0 = mean(X(1,:));
y0 = mean(Y(:,1));

beta = 0.1;
delta = 0.2;
phi = @(x,y) beta*exp(-((X-x0).^2+(Y-y0).^2)/delta);
PHI = phi(X,Y);
for i = 2:1:N+1
    PHI(:,:,i) = RK4(PHI(:,:,end),dt,@rhs_PasScalar);
end
end

%% Ableitungen fuer beliebige Funktionen/Felder F
% from HA01
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
 
k1=rhs(u0          );
k2=rhs(u0+dt/2* k1 );
k3=rhs(u0+dt/2* k2 );
k4=rhs(u0+dt* k3   );

u1 = u0+dt/6 * (k1 + 2*k2 +2*k3+k4);
end