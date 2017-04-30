%MECE 5397 - Computing For Engineers Final Project
clear all; clc;

%Initial Values That Were Given In The Problem Statement That Dictate The
%Bounds For The Problem
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;

%The Number Of Nodes Used In Computing Of The Solution And The
%Incrementation Value That Was Used
n = 64;
m = 64;
dx = 1/n;
dy = 1/m;

%Computes The Coordinate Discretization Values That Are Used In The
%Solution And Places Them Into A 2-D Meshgrid To Be Called Upon Later
x = linspace(ax,bx,n+2);
y = linspace(ay,by,m+2);
[X,Y] = meshgrid(x,y);
F = sin(pi*(X - ax)/(bx - ax)).*cos(pi/2*(2*((Y - ay)/(by - ay)) + 1));

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

U = zeros(n+2,m+2);
U(:,1) = Fb;
U(:,end) = Gb;
U(end,:) = U_b;
ER = 10;

%The Number Of Iterations Used
for Q = 1:1000   
%Guass-Sidel
for k = 2:m+1
    for j = 2:n+1
        Q(j,k) = U(j,k);
        U(j,k) = (-1/C)*(-A*U(j+1,k) - B*U(j,k-1) - A*U(j+1,k) - B*U(j,k+1) - F(j,k));
        relER = abs((U(j,k) - Q(j,k))/U(j,k));
    end
end
%Top-Side Boundary Condition
for j = 2:n+1
        U(1,j) = (-1/C)*(-2*B*U(2,j) - A*U(1,j-1) - A*U(1,j+1) - F(1,j));
end
end
 surf(U)