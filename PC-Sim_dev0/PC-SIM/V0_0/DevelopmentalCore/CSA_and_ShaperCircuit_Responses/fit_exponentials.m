function B = fit_exponentials(I, t, tau)

N = length(tau);

A = zeros(length(t), N);

for i = 1:N
    A(:,i) = exp(-t'/tau(i));
end

B = A \ I(:);

end