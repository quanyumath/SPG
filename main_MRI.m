clear; close all
clc
%--------Ò½Ñ§Í¼Æ¬
[imaVOL,scaninfo] = loadminc('t1_icbm_normal_1mm_pn0_rf20.mnc');
ZZ=imaVOL/max(imaVOL(:));%imshow(ZZ(:,:,38)) 38  88
Z=ZZ(:,:,38);
NN=size(Z);
p=0.9;
Omega = find(rand(prod(NN),1)<p);
X=zeros(size(Z));
X(Omega)=Z(Omega);B=Z(Omega);
%----Ìí¼ÓÔëÉù-------------
B=0.9*imnoise(B,'gaussian',0,0.0001)+0.1*imnoise(B,'gaussian',0,0.01);
% X=zeros(size(Z));X(Omega)=B;imshow(X)
%-------------------
xs=Z;[m,n]=size(Z);sr=p;A=Omega;b=B;
%-------------------
%% Ours(SPG-M)ST
%20;50;3;50
opt.maxiter=200;
tic
[X_Ourst,mo]=SPG(xs,A,b,opt);
time_SPG=toc % imshow(X_Ourst)  imshow(X_Ourst-ZZ(:,:,88)+0.5)
psnr_SPG=PSNR(Z,X_Ourst)




















