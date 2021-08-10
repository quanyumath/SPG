clear; close all 
clc
%------ Setting parameters -------
m=150;n = m;sr=0.8; r=30;
xl =unifrnd(-0.1,0.3,m,r); xr = unifrnd(-0.1,0.3,n,r); xs = xl*xr';
NN=size(xs);Z=xs;
Omega = find(rand(prod(NN),1)<sr);
X=zeros(size(xs));
X(Omega)=xs(Omega);B=xs(Omega);
%----Add noise-------------
c=0.2;
B=(1-c)*imnoise(B,'gaussian',0,0.0001)+c*imnoise(B,'gaussian',0,0.1);
A=Omega;b=B;
%-------------------
%% SPG-M
opt.maxiter=200;
tic
[X_Ourst,~]=SPG_matrix(xs,A,b,opt);
time_SPG=toc; 
REL_SPG=REL(Z,X_Ourst);

