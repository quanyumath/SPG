function y =phi(x,nv)
y_bar=min(1,abs(x)/nv);
y=sum(y_bar);
end