% Please cite: "Q.Yu and X.Zhang. A smoothing 
% proximal gradient algorithm for matrix rank minimization problem"

% Written by: Quan Yu, Tianjin University
% Email: quanyu527@163.com

function [X,iter]=SPG_matrix(Z,Omega,B,opts)

NN=size(Z);maxiter=opts.maxiter;N=min(NN(1),NN(2));
%% Initialize essencial variables
%%%%%%%%%%%Random Matrix Completion%%%%%%%
% -------\mu_0
%mu_bar=opts.mu_bar;rho=1.5;sigma=0.8;gamma_inf=0.1;gamma_sup=5;nv=20;lambda=20;alpha=0.8;
% --------Matrix Size
mu_bar=10;rho=1.5;sigma=0.8;gamma_inf=0.1;gamma_sup=5;nv=20;lambda=20;alpha=0.8;
%%%%%%%%%%%Image Inpainting%%%%%%%%%%%%%
% ------GMM noise
%mu_bar=8;rho=1.5;sigma=0.8;gamma_inf=0.1;gamma_sup=8;nv=20;lambda=20;alpha=0.8;
%----Gaussian noise
%mu_bar=1;rho=1.5;sigma=0.8;gamma_inf=0.1;gamma_sup=8;nv=20;lambda=20;alpha=0.8;
%%%%%%%%%%% MRI Volume Dataset %%%%%%%%%%%%
%mu_bar=1;rho=1.5;sigma=0.8;gamma_inf=0.1;gamma_sup=8;nv=20;lambda=20;alpha=0.8;
X=rand(size(Z));mu=mu_bar;mu_old=mu;mu0=mu_bar;
[~,Sigma,~]=svd(X,'econ');x=diag(Sigma);
fprintf('SPG b\b\b\b Iteration:     ');
for iter=1:maxiter
   fprintf('\b\b\b\b\b%5i',iter); 
   gamma=1/3*gamma_inf+1/3*gamma_sup;
   Xold_Sigma=x;
   %-----Update d(Calculate the number of d==2)------------
   d=sum(x>=nv);
%    d=ones(1,N);d(diag(Sigma)>=nv)=2;
%    for i=1:1:N
%        if Sigma(i,i)<nv
%            d(i)=1;
%        else
%            d(i)=2;
%        end
%    end
   
   %% ------Update X^{k+1}------
   %---------Update gradF-------------
   mb=X(Omega)-B;%Ax-b
   G=gradF(Omega,mb,NN(1),NN(2),mu);
   F_Xold=F(mb,mu_old);
   if mu_old==mu
       F_Xo=F_Xold;
   else
       F_Xo=F(mb,mu);
   end
   %---------Update x^k--------------
   for j=1:1:150
   W=X-mu*G/gamma;
   [U,Sigma,V]=svd(W,'econ');
   w=diag(Sigma);
   tau=lambda*mu/gamma;
   w_bar=w;w_bar(1:d)=w(1:d)+tau/nv;
%    for i=1:1:N
%        if d(i)==1
%            w_bar(i)=w(i);
%        else
%            w_bar(i)=w(i)+tau/nv;
%        end
%    end
%    for i=1:1:N
%        if w_bar(i)<=tau/nv
%            y(i)=0;
%        else 
%            y(i)=w_bar(i)-tau/nv;
%        end
%    end
   x=max(w_bar-tau/nv,0);
   X_hat=U*diag(x)*V';
   DX=X_hat-X;
   F_X=F(X_hat(Omega)-B,mu);
   if F_X<=F_Xo+DX(:)'*G(:)+DX(:)'*DX(:)*gamma/4/mu;
       break;
   end
   gamma=gamma*rho;
   end
   %x(rank_X+1:N)=0; X_hat=U*diag(x)*V';
   X_Sigma=x;X=X_hat;
   if F_X+lambda*phi(X_Sigma,nv)+0.5*mu-F_Xold-lambda*phi(Xold_Sigma,nv)-0.5*mu_old>-alpha*mu^2
        mu_old=mu;mu=mu0/(iter+1)^sigma;
   else
       mu_old=mu;
   end
   if norm(DX(:))/norm(X(:))<1e-3
       break;
   end
   rel_deltaX(iter)=sqrt(norm(X(:)-Z(:))/numel(Z));
end