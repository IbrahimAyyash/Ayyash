%MECE 5397 - Computing For Engineers Final Project
close all; clc;

%Initial Values That Were Given In The Problem Statement That Dictate The
%Bounds For The Problem
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;

%The Number Of Nodes Used In Computing Of The Solution And The
%Incrementation Value That Was Used
n = 30;
dx = 1/n;
dy = 1/n;

%Computes The Coordinate Discretization Values That Are Used In The
%Solution 
x = linspace(ax,bx,n+2);
y = linspace(ay,by,n+2);

%This For Loop Computes The Value Of F And Stores Them In A Matrix Of Size
%n+2 by n+2 That Will Be Called Upon In The Gauss-Sidel and SOR
%Computations
for p = 1:n+2
    for q = 1:n+2
        F(q,p) = sin(pi*(x(q) - ax)/(bx - ax)).*cos(pi/2*(2*((y(p) - ay)/(by - ay)) + 1));
    end
end
%F = 0*F;


%Right-Side Boundary Condition
Gb = ((by - y).^2).*cos((pi*y)/by);
Gby = ((by - ay)^2)*cos((pi*ay)/by);

%Left-Side Boundary Condition
Fb = y.*((by - y).^2);
Fby = ay*((by - ay)^2);

%Bottom-Side Boundary Condition
U_b = Fby + ((x - ax)/(bx - ax))*(Gby - Fby);

%Coefficients Of The Node Values 
A = 1/(dx^2);
B = 1/(dy^2);
C = ((2/dx^2)+(2/dy^2));

%Setting Up A Global Solution Matrix The Contains All The Boundary
%Conditions
U = ones(n+2,n+2);
U(:,1) = Fb;                %Left Boundary Side
U(:,end) = Gb;              %Right Boundary Side
U(1,:) = U_b;               %Bottom Boundary Side
%surf(x,y,U)

%This Step Sets Up Counters To Record The Error Values, Iterations Computed
%And The Time Taken To Finish Computing
tic
Guasscounter = 0;
err = 1;
omega = 1;
save('PoissonEquationAP023')
%This While Loop Allows the Iterative Solver To Keep Computing Until The
%Error Is To The Specified Poiint Where The Error Is Too Small To Consider
while err > 1e-6
D = U;
%Top Side Boundary Conditions Are Computed Here Due To The Neumann
%Conditions Imposed On Them, They Require Their Own Computation Based Of
%The Iterative Solver
 for j = 2:n+1
         U(end,j) = (-1/C)*(-F(end,j) - (2*B)*U(end-1,j) - A*U(end,j-1) - A*U(1,j+1) );
 end
%This Is The Actual Solver, Both Gauss-Sidel And SOR Are Present, However
%When You Plug The Multiplier Value Of w=1 The SOR Behaves Like Gauss-Sidel
for k = 2:n+1
    for j = 2:n+1
        U(j,k) = (-1/C)*(- F(j,k) - A*U(j,k-1) - B*U(j-1,k) - A*U(j,k+1) - B*U(j+1,k));
        U(j,k) = omega*U(j,k) + (1 - omega)*D(j,k);
        Guasscounter = Guasscounter + 1;
    end
end
err = max(max(abs((D-U)./D)));
end
toc

%Plots The Surface And Contour Plots Of The Poisson Solution
surfc(x,y,U)
title('Poisson Surface Plot')
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
figure
contour(x,y,U)
title('Poisson Contour')
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
disp('Computing Iterations:');
disp(Guasscounter);
