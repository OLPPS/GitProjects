function V_out = apply_shaper(V_in, params)

V = V_in;
dt = params.dt;

chain = params.shaper.chain;

for k = 1:length(chain)

    stage = chain{k};

    if strcmp(stage.type,'CR')
        V = CR_filter(V, dt, stage.tau);

    elseif strcmp(stage.type,'RC')
        V = RC_filter(V, dt, stage.tau);

    else
        error('Unknown filter type');
    end
end

V_out = V;

end

%{
function V_out = apply_shaper(V_in, dt, tau_CR, n_RC, tau_RC)

V = V_in;

% CR stage(s)
for k = 1:length(tau_CR)
    V = CR_filter(V, dt, tau_CR(k));
end

% RC stage(s)
for k = 1:n_RC
    V = RC_filter(V, dt, tau_RC(k));
end

V_out = V;

end
%}