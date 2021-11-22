function psnr=PSNR(M,X)
[m,n]=size(M);max=norm(M(:),inf)^2;XDM=X-M;
psnr=10*log10(m*n*max/norm(XDM(:),2)^2);
end