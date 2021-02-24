%---------------------------------------
% CFD1_HA_04_Group BB
%   XY            (xxxxxx)
%   XY            (xxxxxx)
%   Attila Kekesi (xxxxxx)
%
%   MATLAB R2020a
%---------------------------------------

% Riemann-Löser für Wasserwellen in 1-D
% Simulation der 1-D Flachwasserwellen mit Hilfe von Riemannlösungen

%%
clear all
%close all
clc
clf

vid = 1;                    % als Video zu speichern: vid = 1
vid_name = 'CFD1_HA04.avi'; % Name des Videos

%% Anfangswerte, Konstante
g = 9.81;           % Gravitationskonst
CFL = 0.5;          % Courant-Friedrichs-Lewy-Zahl
L = 2*pi;           % x-Laenge
Nx = 50;            % Anzahl der Punkte
Nt = 100;           % Zeitschritte
j = 1:1:Nx;         % Schritte
dx = L/Nx;          % delta_x

x0 = L/2;
b = L/12;
d = 0.5;
h0 = 1;

x = dx*(j-1/2);                 % x-Vektor
h = h0 + d*exp(-(x-x0).^2/b);   % h-Vektor
u = zeros(size(x));             % u-Vektro
m = u.*h;                       % m-Vektro
omega = [h;m];

Omega = omega;      % Matrix fuer Zeitschritte
U = u;              % Matrix fuer Zeitschritte
T = zeros(1,1);     % Matrix fuer Zeit
%% Berechnung
for t = 2:1:Nt    % Schleife fuer Zeitdiskretisierung
    h = Omega(1,:,t-1);
    m = Omega(2,:,t-1);
    hneu = zeros(size(x));
    mneu = zeros(size(x));
    uneu = zeros(size(x));
    a = zeros(size(x));
    hsprung = zeros(size(x));
    msprung = zeros(size(x));
    lambdap = zeros(size(x));
    lambdam = zeros(size(x));
    rm = zeros(size(x));
    hstern = zeros(size(x));
    mstern = zeros(size(x));
    omegastern = zeros(2,length(x));
    f = zeros(2,length(x));
    for n = 1:1:Nx  % Schleife fuer Laengendiskretisierung
        if n < Nx   % wegen Periodizitaet
            hsprung(n) = h(n+1) - h(n);
            msprung(n) = m(n+1) - m(n);
        else
            hsprung(n) = h(1) - h(n);
            msprung(n) = m(1) - m(n);
        end
        lambdap(n) = m(n)/h(n) + (g*h(n))^(1/2);
        lambdam(n) = m(n)/h(n) - (g*h(n))^(1/2);
        a(n) = (g*h(n))^(1/2);
        rm(n) = 1/(2*a(n))*(lambdap(n)*hsprung(n) - 1*msprung(n));
        hstern(n) = h(n) + rm(n)*1;
        mstern(n) = m(n) + rm(n)*lambdam(n);
        omegastern(1,n) = hstern(n);
        omegastern(2,n) = mstern(n);
        f(1,n) = mstern(n);
        f(2,n) = mstern(n)^2/hstern(n) + g/2*hstern(n)^2;
    end
    for n = 1:1:Nx  % Schleife fuer Laengendiskretisierung
        amax = max(a);
        dt = CFL*dx/amax;   % ??? ist es richtig ???
        T(t) = dt;
        if n == 1   % wegen Periodizitaet
            hneu(n) = h(n) - dt/dx*(f(1,n) - f(1,end));
            mneu(n) = m(n) - dt/dx*(f(2,n) - f(2,end));
        else
            hneu(n) = h(n) - dt/dx*(f(1,n) - f(1,n-1));
            mneu(n) = m(n) - dt/dx*(f(2,n) - f(2,n-1));
        end
        uneu(n) = mneu(n)/hneu(n);
    end
    Omega(1,:,t) = hneu;
    Omega(2,:,t) = mneu;
    U(1,:,t) = uneu;
end

%% Plot
% Video
if vid == 1
    Video = VideoWriter(vid_name);
    Video.FrameRate = 13;
    open(Video)
end
% Plot
for t = 1:1:Nt
    subplot(2,1,1)
    plot(x,Omega(1,:,t))
    hold on
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('h - x','FontSize',16,'FontWeight','normal')
    legend({'h'},'location','NorthEast')
    text(5.4,1.879,['t = ',num2str(sum(T(1:t)))],'FontSize',15)
    text(5.05,1.879,['N = ',num2str(t-1)],'FontSize',15)
    xlabel('x')
    ylabel('h')
    xlim([0 x(end)])
    ylim([0 2])
    grid on
    grid minor
    hold off

    subplot(2,1,2)
    plot(x,U(1,:,t))
    hold on
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('u - x','FontSize',16,'FontWeight','normal')
    legend({'u'},'location','NorthEast')
    text(5.4,0.879,['t = ',num2str(sum(T(1:t)))],'FontSize',15)
    text(5.05,0.879,['N = ',num2str(t-1)],'FontSize',15)
    xlabel('x')
    ylabel('u')
    xlim([0 x(end)])
    ylim([-1 1])
    grid on
    grid minor
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