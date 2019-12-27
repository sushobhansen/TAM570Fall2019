clear all;
clc;
format shorte;
n = 50; %number of roots
n1 = 9;
L = 20;
range = [0,1]; %range of r

[roots,iters] = besselj0roots(n);

nquads = 2.^[1:10]';

%Test for qudrature convergence
quadconvergence

%Plot Tr0, Tpr0, and temperature field
Tplots

%Convergence of Tmax
Tmaxconvergence