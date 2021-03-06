%Manufactured solution
close all; clc;
%Initial Values That Were Given In The Problem Statement That Dictate The
%Bounds For The Problem
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;

%The Number Of Nodes Used In Computing Of The Solution And The
%Incrementation Value That Was Used
n = 40;
m = 40;
dx = 1/n;
dy = 1/m;

%Computes The Coordinate Discretization Values That Are Used In The
%Solution And Places Them Into A 2-D Meshgrid To Be Called Upon Later
x = linspace(ax,bx,n+2);
y = linspace(ay,by,m+2);

%These Are The Coefficeints For The Manufactured Solution Function, The
%Simpleest Plot Was Created With All Values Equal To One
a = 1;
b = 1;
w = 1;
save('ManufacturedSol')
%This For Loop Computes The Value Of F And Stores Them In A Matrix Of Size
%n+2 by n+2 That Will Be Called Upon In The Gauss-Sidel and SOR
%Computations
for p = 1:n+2
    for q = 1:n+2
        F(q,p) = ((a^2)+(w^2))*(sin(a*x(p))*cos(w*(y(q)-b)));
    end
end
%Coefficients Of The Node Values
A = 1/(dx^2);
B = 1/(dy^2);
C = ((2/dx^2)+(2/dy^2));

%Setting Up A Global Solution Matrix The Contains All Nodular Solutions
U = ones(n+2,m+2);
U(:,1) = 0;                           %Left Side Boundary
U(:,end) = sin(2*pi*a)*cos(w*(y-b));  %Right Side Boundary
U(1,:) = sin(a*x)*cos(-w*b);          %Bottom Side Boundary

%This Step Sets Up Counters To Record The Error Values, Iterations Computed
%And The Time Taken To Finish Computing
Gausscounter = 0;
err = 1;
omega = 1;

%This Step Is Used In The Top Boundary Condition, It Is An Addition To The
%Ghost Node
z = -w*sin(a*x)*sin(w*((2*pi)-b));

%This While Loop Allows the Iterative Solver To Keep Computing Until The
%Error Is To The Specified Poiint Where The Error Is Too Small To Consider
while err > 1e-8
    D = U;
%This Is The Actual Solver, Both Gauss-Sidel And SOR Are Present, However
%When You Plug The Multiplier Value Of w=1 The SOR Behaves Like Gauss-Sidel
    for k = 2:m+1
        for j = 2:n+1
            U(j,k) = (-1/C)*(- F(j,k) - A*U(j,k-1) - B*U(j-1,k) - A*U(j,k+1) - B*U(j+1,k));
            U(j,k) = omega*U(j,k) + (1 - omega)*D(j,k);
            Gausscounter = Gausscounter + 1;
        end
%Top Side Boundary Conditions Are Computed Here Due To The Neumann
%Conditions Imposed On Them, They Require Their Own Computation Based Of
%The Iterative Solver
            U(end,k) = (-1/C)*(-F(end,k) - B*U(end-1,k) - B*(U(end-1,k)+(2*z(k)*dy)) - A*U(end,k-1) - A*U(end,k+1) );
    end
    err = max(max(abs((D-U)./D)));
end

%Plots The Surface And Contour Plots Of The Poisson Solution
surfc(x,y,U)
title('Discretized Model')
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
disp('Computing Iterations:');
disp(Gausscounter);

%This Step Is Computing The Analytical Soltion and Storing The Values In A
%Matrix Of n+2 by n+2 That Will Be Used To Plot And Compare To DisCrete
%Version
for p = 1:n+2
    for q = 1:n+2
        R(q,p) = sin(a*x(p))*cos(w*(y(q)-b));
    end
end
% R(:,1) = 0;
% R(:,end) = sin(2*pi*a)*cos(w*(y-b));
% R(end,:) = sin(a*x)*cos(w*b);
% R(1,:) = -w*sin(a*x)*sin(w*((2*pi)-b));
%This Step Plots The Analytical Solution In The Form Of A Surface Plot For
%Comparison
figure
surfc(x,y,R)
title('Analytical Model')
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
