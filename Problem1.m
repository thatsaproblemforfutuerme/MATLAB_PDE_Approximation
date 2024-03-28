clc
clear variables
close all

%% 1.e
T =5;%final time
xdom = [0,2];tdom = [0,T];
fi = @(x) 2.5*sin(2*pi*x);
c = .81;
D = .01;
h = .02;k = .001;%steps
Xpoints = xdom(2)/h;Tpoints = tdom(2)/k;
[Xt,Tt, state] = advectionDiffusionEquation1DFTCS(c,D,fi,xdom,tdom,Xpoints,Tpoints,T);
figure(1)
plot(Xt,state(:,end))
grid on
title("1(e)(FTCS) Advection Diffusion at time = "+ T+" seconds" )
xlabel("space")
ylabel("state of fluid")

%% 1.f
%let sigma > .5
k = .1;
Tpoints = tdom(2)/k;%pointless algebra
[Xt,Tt, state] = advectionDiffusionEquation1DFTCS(c,D,fi,xdom,tdom,Xpoints,Tpoints,T);
figure(2)
plot(Xt,state(:,end))
grid on
title("Advection Diffusion at time = "+ T+" seconds when sigma>.5" )
xlabel("space")
ylabel("state of fluid")

%% 2.f 

% k = .01;
% Tpoints = tdom(2)/k;
% [XtDB,TtDB, stateDB] = advectionDiffusionEquation1DBTCS(c,D,fi,xdom,tdom,Xpoints,Tpoints,T);
% figure(3)
% plot(XtDB,stateDB(:,end))
% grid on
% title("2(f)(BTCS) Advection Diffusion at time = "+ T+" seconds" )
% xlabel("space")
% ylabel("state of fluid")

T =5;%final time
xdom = [0,2];tdom = [0,T];
fi = @(x) 2.5*sin(2*pi*x);
c = .81;
D = .01;
h = .02;k = .001;%steps
%none given in problem statement 1.e
Xpoints = xdom(2)/h;Tpoints = tdom(2)/k;%pointless algebra

%t0
[XtBT,TtBT, stateBT] = advectionDiffusionEquation1DBTCS(c,D,fi,xdom,tdom,Xpoints,Tpoints,T);
figure(8)
%t = 0
t=0;
plot(XtBT,stateBT(:,1))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*0)*sin(2*pi*(Xt(:)-c*0));
plot(XtBT,eFzero,'-o')
hold off
grid on
title("(BTCS)Advection Diffusion at time = 0 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")

%t = 1.5
figure(9)
tt =1.5/k;
t = 1.5;
plot(XtBT,stateBT(:,tt))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*t)*sin(2*pi*(Xt(:)-c*t));
plot(XtBT,eFzero,'-o')
hold off
grid on
title("(BTCS)Advection Diffusion at time = 1.5 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")

%t = 3.5
figure(10)
tt =3.5/k;
t = 3.5;
plot(XtBT,stateBT(:,tt))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*t)*sin(2*pi*(Xt(:)-c*t));
plot(XtBT,eFzero,'-o')
hold off
grid on
title("(BTCS)Advection Diffusion at time = 3.5 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")

%t = 4
figure(11)

t = 5;
plot(XtBT,stateBT(:,end))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*t)*sin(2*pi*(Xt(:)-c*t));
plot(XtBT,eFzero,'-o')
hold off
grid on
title("(BTCS)Advection Diffusion at time = 5 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")

%% 3.d
T =5;%final time
xdom = [0,2];tdom = [0,T];
fi = @(x) 2.5*sin(2*pi*x);
c = .81;
D = .01;
h = .02;k = .001;%steps
%none given in problem statement 1.e
Xpoints = xdom(2)/h;Tpoints = tdom(2)/k;%pointless algebra

%t0
[XtCN,TtCN, stateCN] = advectionDiffusionEquation1DCN(c,D,fi,xdom,tdom,Xpoints,Tpoints,T);
figure(4)
%t = 0
t=0;
plot(XtCN,stateCN(:,1))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*0)*sin(2*pi*(Xt(:)-c*0));
plot(XtCN,eFzero,'-o')
hold off
grid on
title("(CN)Advection Diffusion at time = 0 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")
%t = 1.5
figure(5)
tt = 1.5/k;
t = 1.5;
plot(XtCN,stateCN(:,tt))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*t)*sin(2*pi*(Xt(:)-c*t));
plot(XtCN,eFzero,'-o')
hold off
grid on
title("(CN)Advection Diffusion at time = 1.5 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")
%t = 3.5
figure(6)
tt = 3.5/k;
t = 3.5;
plot(XtCN,stateCN(:,tt))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*t)*sin(2*pi*(Xt(:)-c*t));
plot(XtCN,eFzero,'-o')
hold off
grid on
title("(CN)Advection Diffusion at time = 3.5 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")
%t = 5
figure(7)
t = 5;
plot(XtCN,stateCN(:,end))
hold on
eFzero(1,:) = 2.5*exp(-4*(pi^2)*D*t)*sin(2*pi*(Xt(:)-c*t));
plot(XtCN,eFzero,'-o')
hold off
grid on
title("(CN)Advection Diffusion at time = 5 seconds" )
xlabel("space")
ylabel("state of fluid")
legend("approximation","solution")

%% Functions 
function [Xt,Tt, state] = advectionDiffusionEquation1DFTCS(c,D,fi,xdom,tdom,Xpoints,Tpoints,T)
%output descretized spatial and temporal domain corresponding with a final
%time T, and the matrix containing the approximate solutions

%test vars
% clc
% clear variables
% T =5;%final time
% xdom = [0,2];tdom = [0,T];
% 
% f = @(x) 2.5*sin(2*pi*x);
% c = .81;
% D = .01;
% h = .02;k = .001;%steps
h = xdom(2)/Xpoints;
k = tdom(2)/Tpoints;

sigma = (D*k)/h^2;
lamda = (c*k)/h;
%function solves for a new col of approximated points 

%function inputs
%c-advection velocity
%D-diffusion coeff
%fi-function giving intitial state of fluid
%xdom - spital domain
%tdom - time domain
%Xpoints
%Tpoints
%T - final time t




%function outputs
%Xt - discretized spatial domain with end points at a specified T
%Tt - discretized temporal domain corresponding to final time T
%State - matrix of approximate values

%spital and temporal steps 
xstep = h; tstep = k;

%PA empty vector to be filled where each col is a different time and each
%row is a point of fluid in which the state is stored
state = zeros(xdom(2)/xstep,tdom(2)/tstep);


l = length(state(:,1));%length of state
%% FTCS
%build mat A
A = diag((1-2*sigma)*ones(l,1));A = A + diag((sigma+lamda/2)*ones(l-1,1),-1);
A = A + diag((sigma-lamda/2)*ones(l-1,1),1);
%insert last elements of sigma
A(l,2) = (sigma-lamda/2);A(1,l-1) = (sigma+lamda/2);

%initial conditions
x = linspace(0,xdom(2), xdom(2)/xstep);
 for i = 1:l
     state(i,1) = fi(x(i));
 end
t = linspace(0,tdom(2),tdom(2)/tstep);
for j = 2:length(t)
    state(:,j) = A*state(:,j-1);
%     figure(1)
%     plot(x,state(:,j))
end

%discretized spatial domain at T, the last col of the state mat
Xt = x;
%discretized temporal domain at T, 
%I dont understand what this means
Tt = linspace(0,T,Tpoints);



end
function [Xt,Tt, state] = advectionDiffusionEquation1DCN(c,D,fi,xdom,tdom,Xpoints,Tpoints,T)
%output descretized spatial and temporal domain corresponding with a final
%time T, and the matrix containing the approximate solutions

%test vars
% clc
% clear variables
% T =5;%final time
% xdom = [0,2];tdom = [0,T];
% 
% f = @(x) 2.5*sin(2*pi*x);
% c = .81;
% D = .01;
% h = .02;k = .001;%steps
h = xdom(2)/Xpoints;
k = tdom(2)/Tpoints;

sigma = (D*k)/h^2;
lamda = (c*k)/h;
%function solves for a new col of approximated points 

%function inputs
%c-advection velocity
%D-diffusion coeff
%fi-function giving intitial state of fluid
%xdom - spital domain
%tdom - time domain
%Xpoints
%Tpoints
%T - final time t




%function outputs
%Xt - discretized spatial domain with end points at a specified T
%Tt - discretized temporal domain corresponding to final time T
%State - matrix of approximate values

%spital and temporal steps 
xstep = h; tstep = k;

%PA empty vector to be filled where each col is a different time and each
%row is a point of fluid in which the state is stored
state = zeros(xdom(2)/xstep,tdom(2)/tstep);


l = length(state(:,1));%length of state
%% FTCS
%build mat A
A = diag((2+2*sigma)*ones(l,1));A = A + diag((-1*sigma-lamda/2)*ones(l-1,1),-1);
A = A + diag((-1*sigma + lamda/2)*ones(l-1,1),1);
%build mat B
B = diag((2-2*sigma)*ones(l,1));B = B + diag((sigma+lamda/2)*ones(l-1,1),-1);
B = B + diag((sigma-lamda/2)*ones(l-1,1),1);
%insert last elements of sigma
A(l,2) = (-sigma+lamda/2);A(1,l-1) = (-sigma-lamda/2);
B(l,2) = (sigma-lamda/2);B(1,l-1) = (sigma+lamda/2);

%initial conditions
x = linspace(0,xdom(2), xdom(2)/xstep);
 for i = 1:l
     state(i,1) = fi(x(i));
 end
t = linspace(0,tdom(2),tdom(2)/tstep);
for j = 2:length(t)
    state(:,j) = inv(A)*B*state(:,j-1);
%     figure(1)
%     plot(x,state(:,j))
end

%discretized spatial domain at T, the last col of the state mat
Xt = x;
%discretized temporal domain at T, 
%I dont understand what this means
Tt = linspace(0,T,Tpoints);



end
function [Xt,Tt, state] = advectionDiffusionEquation1DBTCS(c,D,fi,xdom,tdom,Xpoints,Tpoints,T)
%output descretized spatial and temporal domain corresponding with a final
%time T, and the matrix containing the approximate solutions

%test vars
% clc
% clear variables
% T =5;%final time
% xdom = [0,2];tdom = [0,T];
% 
% f = @(x) 2.5*sin(2*pi*x);
% c = .81;
% D = .01;
% h = .02;k = .001;%steps
h = xdom(2)/Xpoints;
k = tdom(2)/Tpoints;

sigma = (D*k)/h^2;
lamda = (c*k)/h;
%function solves for a new col of approximated points 

%function inputs
%c-advection velocity
%D-diffusion coeff
%fi-function giving intitial state of fluid
%xdom - spital domain
%tdom - time domain
%Xpoints
%Tpoints
%T - final time t




%function outputs
%Xt - discretized spatial domain with end points at a specified T
%Tt - discretized temporal domain corresponding to final time T
%State - matrix of approximate values

%spital and temporal steps 
xstep = h; tstep = k;

%PA empty vector to be filled where each col is a different time and each
%row is a point of fluid in which the state is stored
state = zeros(xdom(2)/xstep,tdom(2)/tstep);


l = length(state(:,1));%length of state
%% BTCS
%build mat A
A = diag((1+2*sigma)*ones(l,1));A = A + diag((-1*sigma-lamda/2)*ones(l-1,1),-1);
A = A + diag((lamda/2-sigma)*ones(l-1,1),1);
%insert last elements of sigma
A(l,2) = (-1*sigma+lamda/2);A(1,l-1) = (-1*sigma-lamda/2);

%initial conditions
x = linspace(0,xdom(2), xdom(2)/xstep);
 for i = 1:l
     state(i,1) = fi(x(i));
 end
t = linspace(0,tdom(2),tdom(2)/tstep);
for j = 2:length(t)
    state(:,j) = inv(A)*state(:,j-1);
%     figure(1)
%     plot(x,state(:,j))
end

%discretized spatial domain at T, the last col of the state mat
Xt = x;
%discretized temporal domain at T, 
%I dont understand what this means
Tt = linspace(0,T,Tpoints);



end








