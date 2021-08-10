% bb=Ax-b
function z =F(bb,mu)
for i=1:1:length(bb)
    if abs(bb(i))>mu
        y(i)=abs(bb(i));
    else
        y(i)=bb(i)^2/2/mu+mu/2;
    end
end
z=sum(y);
end