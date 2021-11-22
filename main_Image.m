clear; close all
clc
addpath(genpath('picture'));
indimgs = [1:20];
% read data and produce mask
i =9;
id = indimgs(i);
pic_name = [ './picture/',num2str(id),'.tiff'];
I = double(imread(pic_name));
Z = I/max(I(:));  % imshow(Z)
%Z=rgb2gray(Z);
NN=size(Z);
p=0.9;
Omega = find(rand(prod(NN),1)<p);
X=zeros(size(Z));
X(Omega)=Z(Omega);B=Z(Omega);
%----Ìí¼ÓÔëÉù-------------
B=0.9*imnoise(B,'gaussian',0,0.001)+0.1*imnoise(B,'gaussian',0,0.1);
%B=imnoise(B,'gaussian',0,0.0001);
% X(Omega)=B;imshow(X)
%-------------------
xs=Z;[m,n]=size(Z);sr=p;A=Omega;b=B;
%-------------------
%% SPG
%20;50;3;50
opt.maxiter=200;
tic
[X_Ourst,mo]=SPG(xs,A,b,opt);
time_SPG=toc % imshow(X_Ourst)  imshow(X_Ourst-ZZ(:,:,88)+0.5)
psnr_SPG=PSNR(Z,X_Ourst)


