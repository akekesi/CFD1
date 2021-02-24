%---------------------------------------
% CFD1_HA_03_Group BB
%   XY            (xxxxxx)
%   XY            (xxxxxx)
%   Attila Kekesi (xxxxxx)
%
%   MATLAB R2020a
%---------------------------------------

% Burgersgleichung
% Simulation der Burgersgleichung

%%
clear all
%close all
clc
clf

vid = 1;                    % als Video zu speichern: vid = 1
vid_name = 'CFD1_HA03.avi'; % Name des Videos

%% Anfangswerte
L = 2*pi;
Nx = 50;            % Anzahl der Punkte
Nt = 100;           % Anzahl der Zeitschritte
dt = 0.02;          % delta_t
h = L/Nx;           % delta_x
x = (h:h:L)';       % diskritisierter Ort
t = (0:dt:Nt*dt);   % diskritisierte Zeit

u = @(x) 1 + sin(x);    % Anfangsfunktion

global params
params.eps = 0.06;
params.lambda = 0.7;
params.D = D_cen2Ord(Nx,h);
params.D2 = D2_zentral(Nx,h);

% fuer Aufgabe f:
x_ = -5:0.1:5;
y_ = -5:0.1:5;
[X,Y] = meshgrid(x_,y_);
Z = X+1i*Y;
F_ = 1+ Z + Z.^2/2 + Z.^3/6 +Z.^4/24;
F = (real(F_).^2 + imag(F_).^2).^(1/2);

%% Berechnung
u = u(x);           % Anfangswerte
U_explEu = u;       % Anfangswerte in Matrix fuer expl. Euler
U_RK4 = u;          % Anfangswerte in Matrix fuer RK4
EW = zeros(Nx,1);   % Matrix fuer Eigenwerte

for n = 2:1:Nt+1
    [rhs,M] = rhsBurgers(U_explEu(:,n-1),0); % t = 0 Platzhalter
    EW(:,n) = eig(dt*full(M));
    [U_explEu(:,n)] = expl_Eu(U_explEu(:,n-1),dt,M);
    [U_RK4(:,n)] = RK4(U_RK4(:,n-1),dt,t(n),@rhsBurgers);
end

kreis = zeros(1,1);
Nk = 100;
for n = 1:1:Nk
    kreis(n) = (exp(1i*2*pi*(n)/Nk)-1);
end

%% Plot
% Video
if vid == 1
    Video = VideoWriter(vid_name);
    Video.FrameRate = 13;
    open(Video)
end
% Plot
for n = 1:1:Nt+1
    subplot(3,2,[1 2])
    plot(x,U_explEu(:,n),'LineWidth',1.1)
    hold on
    plot(x,U_RK4(:,n),'--','LineWidth',1.3)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Burgers-GL','FontSize',16,'FontWeight','normal')
    legend({'exp. Euler', 'Runge-Kutta 4'},'location','NorthEast')
    text(0.1,-0.5,['t = ',num2str(t(n))],'FontSize',15)
    xlabel('x')
    ylabel('u')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    hold off

    subplot(3,2,3)
    plot(kreis,'r')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Stabilitaetsgrenze expl. Euler','FontSize',16,'FontWeight','normal')
    xlabel('dt*Re(EW)')
    ylabel('dt*Im(EW)')
    xlim([-4  3])
    ylim([-4 4])
    grid on
    grid minor
    daspect([1 1 1])

    subplot(3,2,4)
    for m = 1:1:length(EW(:,n))
        p1 = plot(real(EW(m,n)),imag(EW(m,n)),'kx');
        hold on
    end
    p2 = plot(kreis,'r');
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Spektrum expl. Euler','FontSize',16,'FontWeight','normal')
    legend([p1 p2],{'EW (expl. Euler)','Stabilitaetsgrenze'},'location','NorthEast')
    xlabel('dt*Re(EW)')
    ylabel('dt*Im(EW)')
    xlim([-3  1])
    ylim([-0.5 0.5])
    grid on
    grid minor
    daspect([1 1 1])
    hold off

    subplot(3,2,5)
    contour(X,Y,F,[1,1],'ShowText','on','LineColor','r')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Stabilitaetsgrenze Runge-Kutta 4','FontSize',16,'FontWeight','normal')
    xlabel('dt*Re(EW)')
    ylabel('dt*Im(EW)')
    xlim([-4  3])
    ylim([-4 4])
    grid on
    grid minor
    daspect([1 1 1])

    subplot(3,2,6)
    for m = 1:1:length(EW(:,n))
        p3 = plot(real(EW(m,n)),imag(EW(m,n)),'kx');
        hold on
    end
    contour(X,Y,F,[1,1],'ShowText','on','LineColor','r');
    p4 = plot(NaN,'r');
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Spektrum Runge-Kutta 4','FontSize',16,'FontWeight','normal')
    legend([p3 p4],{'EW (Runge-Kutta 4)','Stabilitaetsgrenze'},'location','NorthEast')
    xlabel('dt*Re(EW)')
    ylabel('dt*Im(EW)')
    xlim([-3  1])
    ylim([-0.5 0.5])
    grid on
    grid minor
    daspect([1 1 1])
    hold off
    drawnow
    
    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

if vid == 1
    close(Video)
end

%% Ableitungsmatrix (zentrale Diff. 2.Ord.)
function D = D_cen2Ord(N,h)
e = ones(N,1);
p = [-N+1 -1 1 N-1];
D = 1/(2*h)*spdiags([e -e e -e],p,N,N);
end

%% 2. Abl. (zentrale Diff.)
function D = D2_zentral(N,h)
e = ones(N,1);
p = [-N+1 -1 0 1 N-1];
D = 1/(h^2)*spdiags([e e -2*e e e],p,N,N);
end

%% Rechte Seite der Burgers-GL
function [rhs, M] = rhsBurgers(u,t)
global params;
U = diag(u);
M = -(U*params.D - params.eps*params.D2);
rhs = M*u;
end

%% expl. Euler
function [u] = expl_Eu(u0,dt,M)
u = u0+dt*M*u0;
end

%% Runge-Kutta 4. Ordnung
function [u] = RK4(u0,dt,t0,rhs)
k1 = rhs(u0,t0);
k2 = rhs(u0+k1*dt/2,t0+dt/2);
k3 = rhs(u0+k2*dt/2,t0+dt/2);
k4 = rhs(u0+k3*dt/2,t0+dt);
u = u0 + dt/6*(k1+2*k2+2*k3+k4);
end

%% Antworte
% Zur Aufgabe c): Das Verfahren ist stabil, weil die Eigenwerte von deltat*M
% bei jedem Zeitschritt im betrachteten Zeitintervall innerhalb des 
% Stabilitaetskreises vom expliziten Euler-Verfahren liegen.