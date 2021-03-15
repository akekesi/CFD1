%---------------------------------------
% CFD1_HA_04_Group BB
%   XY            (xxxxxx)
%   XY            (xxxxxx)
%   Attila Kekesi (xxxxxx)
%
%   MATLAB R2020a
%---------------------------------------

% Numerische Ableitung und Übertragunsverhalten
% Numerische Ableitung
% Übertragungsverhalten der diskreten Ableitungen
% Verbessertes Übertragungsverhalten

%%
clear all
%close all
clc
clf

%% Anfangswerte

L = 2*pi;
N = 100;            % gerade
h = L/N;
x = (h:h:L)';
n = 5;              % n = 0,1,2,...,N/2
kh = (0:1:N/2)*h;   % k*h = 0,...,pi
kn = (2*pi*n)/L;
a3 = -0.2:0.1:0.2;

%% Berechnung

% Aufgabe_1
y = sin(x);
y2 = [y;y];
D = D_cen2Ord(N,h);
dy_ana = cos(x);
dy_num = D*y;

% Aufgabe_2
u = exp(1i*kn*x);
du = D*u;
kn_ = du./(1i*u);   % modifizierte Wellenzahl zu kn
Kn_ = kMod(D,L);    % modifizierte Wellenzahl zu alle kn

% Aufgabe_3
Kn_2 = zeros(size(Kn_));
for m = 1:1:length(a3)
    D2 = D_central(N,h,a3(m));
    Kn_2(m,:) = kMod(D2,L);
end

DTW = D_central(N,h,0.15912);
Kn_TW = kMod(DTW,L);

A = [ 1 1 1; 1 4 9;1 16 81];
b = [1 0 0];
a = A\b';
D6 = D_central(N,h,a(3));
Kn_6 = kMod(DTW,L);

%% Plot
% Aufgabe_1
subplot(2,3,1)
p1 = plot(x,y);
hold on
p2 = plot(x,dy_ana);
p3 = plot(x,dy_num,'kx');
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('x - sin(x)','FontSize',16,'FontWeight','normal')
legend([p1 p2 p3],["y","y' (ana)","y' (num)"],'location','NorthEast')
xlabel('x')
ylabel('sin(x)')
xlim([0 L])
xticks([0 1/2*pi pi 3/2*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
grid on
grid minor
hold off

subplot(2,3,4)
plot(y2)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('j - sin(x)','FontSize',16,'FontWeight','normal')
legend('sin(x)','location','NorthEast')
xlabel('j')
ylabel('sin(x)')
xlim([0 2*length(y)-1])
grid on
grid minor

% Aufgabe_2
subplot(2,3,2)
p1 = plot(real(kn_));
hold on
p2 = plot(imag(kn_));
hp = fix(2*real(kn_(1)));
plot(hp)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title({
    'j - Re(k''_n),  j - Im(k''_n)'
    ['n = ',num2str(n)]
    },'FontSize',16,'FontWeight','normal')
legend([p1 p2],{'Re(k''_n)','Im(k''_n)'},'location','NorthEast')
xlabel('j')
ylabel('Re(k''_n),  Im(k''_n)')
yticks(-1:1:hp)
grid on
grid minor
hold off

subplot(2,3,5)
p1 = plot(kh,kh);
hold on
p2 = plot(kh,real(Kn_)*h);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('kh - kh,  kh - Re(k''_n)h','FontSize',16,'FontWeight','normal')
legend([p1 p2],{'kh','Re(k''_n)h'},'location','NorthEast')
xlabel('kh')
ylabel('kh,  Re(k''_n)h')
xticks([0 1/2*pi pi])
xticklabels({'0','\pi/2','\pi'})
yticks([0 1/2*pi pi])
yticklabels({'0','\pi/2','\pi'})
grid on
grid minor
hold off

% Aufgabe_3
subplot(2,3,3)
p1 = plot(kh,kh);
hold on
p2 = plot(kh,real(Kn_TW)*h);
p3 = plot(kh,real(Kn_6)*h);
legend([p1 p2 p3],{'kh','TW','6.'},'location','NorthEast')
for m = 1:1:length(a3)
    legend_a3 = ['\alpha_3 = ',num2str(a3(m))];
    plot(kh,real(Kn_2(m,:))*h,'DisplayName',legend_a3);
    hold on
end
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('kh - kh,  kh - Re(k''_n(\alpha_3))h','FontSize',16,'FontWeight','normal')
xlabel('kh')
ylabel('kh,  Re(k''_n(\alpha_3))h')
xticks([0 1/2*pi pi])
xticklabels({'0','\pi/2','\pi'})
yticks([0 1/2*pi pi])
yticklabels({'0','\pi/2','\pi'})
grid on
grid minor
hold off

txt = "Tam und Webb hat das beste Übertragungsverhalten";
annotation('textbox', [0.7, 0.3, 0.1, 0.1], 'String', txt, 'FitBoxToText','on')


%% Ableitungsmatrix (zentrale Diff. 2.Ord.)
function D = D_cen2Ord(N,h)
e = ones(N,1);
p = [-N+1 -1 1 N-1];
D = 1/(2*h)*spdiags([e -e e -e],p,N,N);
end

%% Ableitungsmatrix (zentrale Diff. 4.Ord.)
function D = D_central(N,h,a3)
A = [1 1; 1  4];
b = [1-a3 -9*a3];
a12 = A\b';
a1 = a12(1)/(2*h);
a2 = a12(2)/(4*h);
a3 = a3/(6*h);
e = ones(N,1);
p = [-N+3 -N+2 -N+1 -3 -2 -1 1 2 3 N-1 N-2 N-3];
D = spdiags([a3*e a2*e a1*e -a3*e -a2*e -a1*e a1*e a2*e a3*e -a1*e -a2*e -a3*e],p,N,N);
end

%% kMod(D,L)
function Kn_ = kMod(D,L)
N = length(D);
h = L/N;
x = (h:h:L)';
Kn_ = zeros(1,1);
for n = 0:1:N/2
    kn = (2*pi*n)/L;
    u = exp(1i*kn*x);
    du = D*u;
    Kn_(n+1) = du(1)/(1i*u(1));
end
end