function rel=REL(Z,X)
rel=norm(X(:)-Z(:))/sqrt(numel(Z));
end
