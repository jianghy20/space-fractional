function result = coef_matrix(N,gam,c,h)
result = zeros(N-1);
for j = 1:N-1
    for k = 1:N-1
        result(j,k) = c(j-k,gam)/(h^gam);
    end
end
