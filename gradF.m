    function y =gradF(AA,bb,m,n,mu)
        % y = A'*b
        y = zeros(m*n,1);
        for i=1:1:length(bb)
            if bb(i)>=mu
                cb(i)=1;
            elseif bb(i)<=-mu
                cb(i)=-1;
            else
                cb(i)=bb(i)/mu;
            end
        end
        y(AA) = cb;
        y = reshape(y,m,n);
    end

