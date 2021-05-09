%---------------------------------------
% CFD1_HA_04_Group BB
%   XY            (xxxxxx)
%   XY            (xxxxxx)
%   Attila Kekesi (xxxxxx)
%
%   MATLAB R2020a
%---------------------------------------

% Einfache Zeitintegrationen
% Simulation der Transportgleichung
% Verschiedene Einflüsse der Ableitung

%%
clear all
%close all
clc
clf

vid = 1;                    % als Video zu speichern: vid = 1
vid_name = 'CFD1_HA02.avi'; % Name des Videos

%% Anfangswerte, Konstante
lambda = 1;
epsilon_1 = 0.01;
epsilon_2 = 0;
a3_1 = 0;
a3_2 = 0.15912;
L = 2*pi;
Nx = 100;
h = L/Nx;
x = (h:h:L)';
dt = h/lambda;
Nt = 100;
f1 = @(x) 1+sin(x);
f2 = @(N) [(1:N/2)/(N/2) 1-(1:N/2)/(N/2)];

%% Berechnung
D = D_central(Nx,h,a3_1);
D2 = D2_zentral(Nx,h);
DL = D_links(Nx,h);
DR = D_rechts(Nx,h);
M_1Z = M_matrix(lambda,D,epsilon_1,D2);
M_2L = M_matrix(lambda,DL,epsilon_1,D2);
M_2R = M_matrix(lambda,DR,epsilon_1,D2);
%full(D)
%full(D2)
%full(DL)
%full(DR)
%full(M_1)
EW_1Z = eig(dt*full(M_1Z));
EW_2L = eig(dt*full(M_2L));
EW_2R = eig(dt*full(M_2R));

% 1: expl. Euler 1.Abl. zentral 4.Ord. 
% 2: impl. Euler 1.Abl. zentral 4.Ord.
% 3: Trapezverf. 1.Abl. zentral 4.Ord.
% 4: expl. Euler 1.Abl. links 2.Ord.
% 5: impl. Euler 1.Abl. links 2.Ord.
% 6: Trapezverf. 1.Abl. links 2.Ord.
% 7: expl. Euler 1.Abl. rechts 2.Ord.
% 8: impl. Euler 1.Abl. rechts 2.Ord.
% 9: Trapezverf. 1.Abl. rechts 2.Ord.
% 21: Trapezverf. 1.Abl. zentral 2.Ord.
% 22: Trapezverf. 1.Abl. Tam Webb


F_expl_Eu_1 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,1);
F_expl_Eu_2 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,1);
F_impl_Eu_1 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,2);
F_impl_Eu_2 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,2);
F_Trapez_1 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,3);
F_Trapez_2 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,3);

F_expl_Eu_3 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,4);
F_expl_Eu_4 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,4);
F_impl_Eu_3 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,5);
F_impl_Eu_4 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,5);
F_Trapez_3 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,6);
F_Trapez_4 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,6);

F_expl_Eu_5 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,7);
F_expl_Eu_6 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,7);
F_impl_Eu_5 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,8);
F_impl_Eu_6 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,8);
F_Trapez_5 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_1,9);
F_Trapez_6 = Berechnung(f1,x,Nx,h,dt,Nt,lambda,epsilon_2,a3_1,9);

F_Trapez_7 = Berechnung(f2,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_2,21);
F_Trapez_8 = Berechnung(f2,x,Nx,h,dt,Nt,lambda,epsilon_1,a3_2,22);
%% Antwort
% 1/a:
%   Nach ein paar Zeitschritten wird der Berechnunginstabil
%   Wann stabil?
% 1/b,c:
%   Fuer epsilon = 0 gibt es keine Reibung (ausser num. Reibung)

%% Plot
% Video
if vid == 1
    Video = VideoWriter(vid_name);
    Video.FrameRate = 13;
    open(Video)
end
% Plot
subplot(3,5,[1 6 11])
for n = 1:1:length(EW_1Z)
    p1 = plot(real(EW_1Z(n)),imag(EW_1Z(n)),'kx');
    hold on
    p2 = plot(real(EW_2L(n)),imag(EW_2L(n)),'rx');
    p3 = plot(real(EW_2R(n)),imag(EW_2R(n)),'bx');
end
kreis = zeros(1,1);
Nk = 100;
for n = 1:1:Nk
    kreis(n) = (exp(1i*2*pi*(n)/Nk)-1);
end
p4 = plot(kreis,'k');
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('Spektrum','FontSize',16,'FontWeight','normal')
legend([p1 p2 p3 p4],{'zentral','links','rechts','Kreis'},'location','NorthEast')
xlabel('dt*Re(EW)')
ylabel('dt*Im(EW)')
grid on
grid minor
daspect([1 1 1])
if vid == 1
    frame = getframe(gcf);
    writeVideo(Video,frame);
end

subplot(3,5,2)
for n = 1:1:Nt+1
    plot(x,F_expl_Eu_1(n,:))
    hold on
    plot(x,F_expl_Eu_2(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('expl. Euler (zentral)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,7)
for n = 1:1:Nt+1
    plot(x,F_impl_Eu_1(n,:))
    hold on
    plot(x,F_impl_Eu_2(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('impl. Euler (zentral)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,12)
for n = 1:1:Nt+1
    plot(x,F_Trapez_1(n,:))
    hold on
    plot(x,F_Trapez_2(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Trapez (zentral)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,3)
for n = 1:1:Nt+1
    plot(x,F_expl_Eu_3(n,:))
    hold on
    plot(x,F_expl_Eu_4(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('expl. Euler (links)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,8)
for n = 1:1:Nt+1
    plot(x,F_impl_Eu_3(n,:))
    hold on
    plot(x,F_impl_Eu_4(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('impl. Euler (links)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,13)
for n = 1:1:Nt+1
    plot(x,F_Trapez_3(n,:))
    hold on
    plot(x,F_Trapez_4(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Trapez (links)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,4)
for n = 1:1:Nt+1
    plot(x,F_expl_Eu_5(n,:))
    hold on
    plot(x,F_expl_Eu_6(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('expl. Euler (rechts)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,9)
for n = 1:1:Nt+1
    plot(x,F_impl_Eu_5(n,:))
    hold on
    plot(x,F_impl_Eu_6(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('impl. Euler (rechts)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,14)
for n = 1:1:Nt+1
    plot(x,F_Trapez_5(n,:))
    hold on
    plot(x,F_Trapez_6(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('Trapez (rechts)','FontSize',16,'FontWeight','normal')
    legend({['esp = ',num2str(epsilon_1)],['esp = ',num2str(epsilon_2)]},'location','NorthEast')
    xlabel('x')
    ylabel('1+sin(x)')
    xlim([0 x(end)])
    ylim([-1 3])
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

subplot(3,5,[5 10 15])
for n = 1:1:Nt+1
    plot(F_Trapez_7(n,:))
    hold on
    plot(F_Trapez_8(n,:),"--")
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title('zentral 2.Ord. / Tam Webb','FontSize',16,'FontWeight','normal')
    legend([" zentral 2.Ord.","Tam Webb"],'location','NorthEast')
    xlabel('x')
    ylabel('y')
    xlim([0 Nx])
    ylim([-0.5 1.5])
    grid on
    grid minor
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

if vid == 1
    close(Video)
end

%% 2. Abl. (zentrale Diff.)
function D = D2_zentral(N,h)
e = ones(N,1);
p = [-N+1 -1 0 1 N-1];
D = 1/(h^2)*spdiags([e e -2*e e e],p,N,N);
end

%% 1. Abl. (zentrale Diff. 4.Ord.)
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

%% Ableitungsmatrix (zentrale Diff. 2.Ord.)
function D = D_cen2Ord(N,h)
e = ones(N,1);
p = [-N+1 -1 1 N-1];
D = 1/(2*h)*spdiags([e -e e -e],p,N,N);
end

%% 1. Abl. (linke Diff.)
function D = D_links(N,h)
e = ones(N,1);
p = [-1 0 N-1];
D = 1/(h)*spdiags([-e e -e],p,N,N);
end

%% 1. Abl. (rechte Diff.)
function D = D_rechts(N,h)
e = ones(N,1);
p = [-N+1 0 1];
D = 1/(h)*spdiags([e -e e],p,N,N);
end

%% M = -lambda*D + epsilon*D2;
function M = M_matrix(lambda,D,epsilon,D2)
M = -lambda*D + epsilon*D2;
end
%% expl. Euler
function f = expl_Eu(f,dt,M)
f = f.'+dt*M*f.';
end

%% impl. Euler
function f = impl_Eu(f,dt,M)
f = (eye(size(M)) - dt*M)\f.';
end

%% Trapezverf.
function f = Trapez(f,dt,M)
f = (eye(size(M)) - dt/2*M)\((eye(size(M)) + dt/2*M)*f.');
end

%% Berechnung
function F = Berechnung(f,x,Nx,h,dt,Nt,lambda,epsilon,a3,Wahl)
if Wahl == 1 || Wahl == 2 || Wahl == 3
    D = D_central(Nx,h,a3);
    D2 = D2_zentral(Nx,h);
    M = M_matrix(lambda,D,epsilon,D2);
    F = f(x).';
elseif Wahl == 4 || Wahl == 5 || Wahl == 6
    D = D_links(Nx,h);
    D2 = D2_zentral(Nx,h);
    M = M_matrix(lambda,D,epsilon,D2);
    F = f(x).';
elseif Wahl == 7 || Wahl == 8 || Wahl == 9
    D = D_rechts(Nx,h);
    D2 = D2_zentral(Nx,h);
    M = M_matrix(lambda,D,epsilon,D2);
    F = f(x).';
elseif Wahl == 21
    D = D_cen2Ord(Nx,h);
    D2 = D2_zentral(Nx,h);
    M = M_matrix(lambda,D,epsilon,D2);
    F = f(Nx);
elseif Wahl == 22
    D = D_central(Nx,h,a3);
    D2 = D2_zentral(Nx,h);
    M = M_matrix(lambda,D,epsilon,D2);
    F = f(Nx);
end

for n = 2:1:Nt+1
    if Wahl == 1
        F(n,:) = (expl_Eu(F(n-1,:),dt,M))';
    elseif Wahl == 2
        F(n,:) = (impl_Eu(F(n-1,:),dt,M))';
    elseif Wahl == 3
        F(n,:) = (Trapez(F(n-1,:),dt,M))';
    elseif Wahl == 4
        F(n,:) = (expl_Eu(F(n-1,:),dt,M))';
    elseif Wahl == 5
        F(n,:) = (impl_Eu(F(n-1,:),dt,M))';
    elseif Wahl == 6
        F(n,:) = (Trapez(F(n-1,:),dt,M))';
    elseif Wahl == 7
        F(n,:) = (expl_Eu(F(n-1,:),dt,M))';
    elseif Wahl == 8
        F(n,:) = (impl_Eu(F(n-1,:),dt,M))';
    elseif Wahl == 9
        F(n,:) = (Trapez(F(n-1,:),dt,M))';
    elseif Wahl == 21 || Wahl == 22
        F(n,:) = (Trapez(F(n-1,:),dt,M))';
    end
end
end