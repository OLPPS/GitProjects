function V_csa = CSA_model(I, t, dt, Cf, Rf, mode, I_reset, trace_length)
%trace_length is the desired legnth of the output trace in ns

trace_length = trace_length / 10^9;

nt_in = length(t);
nt_out = round(trace_length/dt);

% ---Extend time ---
t_ext = (0:nt_out-1) * dt;

% ---Extend input current ---
I_ext = zeros(1,nt_out);
I_ext(1:min(nt_in, nt_out)) = I(1:min(nt_in,nt_out));

V_csa = zeros(1,nt_out);

if strcmp(mode,'exp')

    tau_f = Rf*Cf;

    for n=2:nt_out
        V_csa(n) = V_csa(n-1) ...
                 + (dt/Cf)*I_ext(n) ...
                 - (dt/tau_f)*V_csa(n-1);
    end

elseif strcmp(mode,'tri')

    for n=2:nt_out
        V_csa(n) = V_csa(n-1) ...
                 + (dt/Cf)*(I_ext(n) - I_reset);

        if V_csa(n) < 0
            V_csa(n) = 0;   % clamp
        end
    end

else
    error('Unknown CSA mode');
end

end