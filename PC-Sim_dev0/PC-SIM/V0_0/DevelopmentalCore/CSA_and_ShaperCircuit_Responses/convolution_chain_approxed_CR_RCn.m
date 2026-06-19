function V = convolution_chain(Q, t, Rvals, Cvals)

dt = t(2)-t(1);
nt = length(t);

% --- current ---
I = gradient(Q, dt);
I = I(:)';   % force row vector

% --- zero pad to full time ---
I_full = zeros(1, nt);
I_full(1:length(I)) = I;

% --- build impulse response ---
tau_f = Rvals(1)*Cvals(1);
h_total = exp(-t/tau_f);

for k=2:length(Rvals)
    tau_s = Rvals(k)*Cvals(k);
    h_rc = (t./tau_s).*exp(-t/tau_s);
    h_total = conv(h_total, h_rc, 'same');
end

h_total = h_total(:)';

% --- FFT convolution ---
V = real(ifft(fft(I_full) .* fft(h_total))) * dt;

end