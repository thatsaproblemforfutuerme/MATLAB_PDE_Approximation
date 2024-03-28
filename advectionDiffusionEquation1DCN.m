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