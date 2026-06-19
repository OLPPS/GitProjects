function V_out = RC_filter(V_in, dt, tau)

nt = length(V_in);
beta = exp(-dt/tau);

V_out = zeros(1,nt);

for n = 2:nt
    V_out(n) = beta * V_out(n-1) ...
             + (1 - beta) * V_in(n);
end

end