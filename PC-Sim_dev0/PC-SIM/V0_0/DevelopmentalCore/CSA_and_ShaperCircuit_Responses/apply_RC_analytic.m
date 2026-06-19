function V_out = apply_RC_analytic(V_in, t, tau)

dt = t(2)-t(1);

h = exp(-t/tau);

V_out = real(ifft(fft(V_in) .* fft(h))) * dt;

end