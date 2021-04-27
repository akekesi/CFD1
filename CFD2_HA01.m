%---------------------------------------
% CFD2_HA_01_Group L
%   XY            (xxxxxx)
%   XY            (xxxxxx)
%   Attila Kekesi (xxxxxx)
%
%   MATLAB R2020a
%---------------------------------------

% Gitterstreckung: Lokale Basis und Geometriefaktoren
% Einfache 2-dimensionale Gitter

%%
clear all
%close all
clc
clf

%% Anfangswerte, Konstante
global params
alpha = 0.3;
N_xi = 30;
N_eta = 25;
xi_end = 2*pi;
eta_end = 2*pi;
h_xi = xi_end/(N_xi-1);     % N_xi diskrete Werte mit 0?
h_eta = eta_end/(N_eta-1);  % N_eta diskrete Werte mit 0?
xi = (0:h_xi:xi_end);       % xi = linspace(0,2*pi,N_xi);
eta = (0:h_eta:eta_end);    % eta = linspace(0,2*pi,N_eta);
[XI,ETA] = meshgrid(xi,eta);

x1 = @(xi,eta) xi + alpha * sin(xi + eta);
y1 = @(xi,eta) eta + alpha * sin(xi + eta);
X1 = x1(XI,ETA);
Y1 = y1(XI,ETA);

x2 = @(xi,eta) xi + alpha * sin(xi - eta);
y2 = @(xi,eta) eta + alpha * sin(xi - eta);
X2 = x2(XI,ETA);
Y2 = y2(XI,ETA);

Z = zeros(size(XI));

%% Berechnung
D_offen_xi = D_offen(N_xi,h_xi);
D_offen_eta = D_offen(N_eta,h_eta);
params.D_kron_xi = kron(speye(N_eta),D_offen_xi);
params.D_kron_eta = kron(D_offen_eta,speye(N_xi));

% mit X1, Y1
ex1_xi = D_xi(X1);
ey1_xi = D_xi(Y1);
ex1_eta = D_eta(X1);
ey1_eta = D_eta(Y1);
g1_xixi = ex1_xi.^2 + ey1_xi.^2;
g1_etaeta = ex1_eta.^2 + ey1_eta.^2;
g1_xieta = ex1_xi.*ex1_eta + ey1_xi.*ey1_eta;
g1_xieta_norm = g1_xieta.*(g1_xixi.*g1_etaeta).^(-1/2);
J1 = ex1_xi.*ey1_eta - ex1_eta.*ey1_xi;

% mit X2, Y2
ex2_xi = D_xi(X2);
ey2_xi = D_xi(Y2);
ex2_eta = D_eta(X2);
ey2_eta = D_eta(Y2);
g2_xixi = ex2_xi.^2 + ey2_xi.^2;
g2_etaeta = ex2_eta.^2 + ey2_eta.^2;
g2_xieta = ex2_xi.*ex2_eta + ey2_xi.*ey2_eta;
g2_xieta_norm = g2_xieta.*(g2_xixi.*g2_etaeta).^(-1/2);
J2 = ex2_xi.*ey2_eta - ex2_eta.*ey2_xi;

% f = sin(x) mit X1, Y1
f = sin(X1);
dx_f_ana = cos(X1);
dy_f_ana = zeros(size(f));
dx_f_num = (D_xi(f).*ey1_eta - D_eta(f).*ey1_xi)./J1;
dy_f_num = (-D_xi(f).*ex1_eta + D_eta(f).*ex1_xi)./J1;

%% Plot
figure(1)
subplot(2,2,1)
pcolor(X1,Y1,Z)
hold on
p1 = quiver(X1,Y1,ex1_xi,ey1_xi);
p2 = quiver(X1,Y1,ex1_eta,ey1_eta);
title({'Lokale Basis','Gitter-1'},'FontSize',16,'FontWeight','normal')
legend([p1,p2],["e_{\xi}","e_{\eta}"],'location','NorthEast')
xlabel('x')
ylabel('y')
daspect([1 1 1])

subplot(2,2,3)
pcolor(X2,Y2,Z)
hold on
p1 = quiver(X2,Y2,ex2_xi,ey2_xi);
p2 = quiver(X2,Y2,ex2_eta,ey2_eta);
title({'Lokale Basis','Gitter-2'},'FontSize',16,'FontWeight','normal')
legend([p1,p2],["e_{\xi}","e_{\eta}"],'location','NorthEast')
xlabel('x')
ylabel('y')
daspect([1 1 1])

subplot(2,2,2)
pcolor(X1,Y1,g1_xieta_norm)
title({'Nichtdiagonalelemente des Metriktensors','Gitter-1'},'FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
caxis([-1 1])
daspect([1 1 1])

subplot(2,2,4)
pcolor(X2,Y2,g2_xieta_norm)
title({'Nichtdiagonalelemente des Metriktensors','Gitter-2'},'FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
caxis([-1 1])
daspect([1 1 1])

figure(2)
subplot(2,4,[1,5])
pcolor(X1,Y1,f)
title('f = sin(x)','FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
caxis([-1 1])
daspect([1 1 1])

subplot(2,4,2)
pcolor(X1,Y1,dx_f_num)
title('df/dx num.','FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
caxis([-1 1])
daspect([1 1 1])

subplot(2,4,6)
pcolor(X1,Y1,dy_f_num)
title('df/dy num.','FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
%caxis([-0.01 0.01])
daspect([1 1 1])

subplot(2,4,3)
pcolor(X1,Y1,dx_f_ana)
title('df/dx ana.','FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
caxis([-1 1])
daspect([1 1 1])

subplot(2,4,7)
pcolor(X1,Y1,dy_f_ana)
title('df/dy ana.','FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
caxis([-0.01 0.01])
daspect([1 1 1])

subplot(2,4,4)
pcolor(X1,Y1,dx_f_num-dx_f_ana)
title('df/dx num. - df/dx ana.','FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
%caxis([-1 1])
daspect([1 1 1])

subplot(2,4,8)
pcolor(X1,Y1,dy_f_num-dy_f_ana)
title('df/dy num. - df/dx ana.','FontSize',16,'FontWeight','normal')
xlabel('x')
ylabel('y')
colorbar
%caxis([-1 1])
daspect([1 1 1])

%% nicht-periodische zentrale Differenzenmatrix
function D = D_offen(N,h)
e = ones(N,1);
p = [-1 1];
D = spdiags([-e e],p,N,N);
D(1,1:3) = [-3 4 -1];
D(end,end-2:end) = [1 -4 3];
D = 1/(2*h)*D;
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
F = reshape(F',n*m,1);
F_eta = D_kron_eta*F;
F_eta = reshape(F_eta,m,n)';
end