function V_out = CR_filter(V_in, dt, tau)

nt = length(V_in);
alpha = exp(-dt/tau);

V_out = zeros(1,nt);

for n = 2:nt
    V_out(n) = alpha * V_out(n-1) ...
             + V_in(n) ...
             - V_in(n-1);
end

end