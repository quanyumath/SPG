function rel=REL(Z,X)
rel=sqrt(norm(X(:)-Z(:))/numel(Z));
end