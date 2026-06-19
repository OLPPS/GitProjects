%MC Data formatting and set standard file names

%Choose if doing CIE, basic charge cloud projection or lookup table for
%responses (is this different from CIE? perhaps need LUT options? Perhaps use something to generate a LUT?)

%CIE Map formatting and set standard file names

%Generate formatting info and selection info: perhaps including a 
%walkthrough GUIDE Script, like the scintillator type, for PC-SIM to
%read in and format itself with.

%Then--- techinically PC-SIM starts here?
    %INITIALISATION OF SCANNER BUILD
        %Import detector
        %Read in time formatting information
        %Splice MC data to extract relevant timeslice data (looping over all time slices)
     %SIGNAL GENERATION STEPS
            %Perform pixelisation
                %Calculate pixel numbers ffrom global coordinates
                %Convert from global to intrapixel coordinates
            %Perform dimensional scaling (if option needed)
                %Normalise intrapixel coordinates
            %Generate PseudoEvents (CIE or other method, if option needed)
                %Identify sharing events based on boundary rule or LUT if needed (a purely CSE LUT, calculated based on thresholding perhaps, would be useful here).
                %Expand list of events to make space for shared events in ordered sequence
            %Adjust times to account for finite drift time (if option on).
                %Adjust times based on time to signal peaking in isolation (or reaching a set threshold, say 90%?). To account for negative signals, this may need to be a modulus function, i.e. "Regardless of the sign, when is the impact of this siganl on the induced event greatest?" 
            %Perform LUT transformations to generate signal trains or other abstracted signal data needed.
            %Signal train may be stored as red or may be abstracted as several key numbers for a fitted function. (may involve projection ones too)

    %SIGNAL PROCESSING





    "Moved sinograms folder into main section and added template script for the overall flow with PC-SIM and this folder plan."




    %%
    function Vmap = process_detector_signal(Qmap, t, Cf, Rf, Rc, Cc)
% PROCESS_DETECTOR_SIGNAL
% Takes Prettyman Q(t,x,y,z) and applies CSA + CR-RC shaping
%
% Inputs:
%   Qmap (nx, ny, nz, nt) : induced charge
%   t    (1 x nt)         : time vector (seconds)
%   Cf                     : feedback capacitance (F)
%   Rf                     : feedback resistance (Ohm)
%   Rc, Cc                 : shaping CR-RC components
%
% Output:
%   Vmap (nx, ny, nz, nt) : output voltage signal

[nx, ny, nz, nt] = size(Qmap);

dt = t(2) - t(1);

% --- Step 1: Compute current I = dQ/dt ---
Imap = diff(Qmap,1,4) / dt;

% Pad to keep same size
Imap(:,:,:,nt) = Imap(:,:,:,nt-1);

% --- Step 2: CSA impulse response ---
tau_f = Rf * Cf;
h_csa = (1/Cf) * exp(-t / tau_f);

% --- Step 3: Shaper impulse response (CR-RC) ---
tau_s = Rc * Cc;

% CR-RC impulse response (approx)
h_shaper = (t./tau_s) .* exp(-t/tau_s);

% Normalize (optional but helps stability)
h_shaper = h_shaper / max(h_shaper);

% --- Step 4: Combined response ---
h_total = conv(h_csa, h_shaper);
h_total = h_total(1:nt); % truncate to original time length

% --- Step 5: Convolution for each voxel ---
Vmap = zeros(nx, ny, nz, nt);

for ix = 1:nx
    for iy = 1:ny
        for iz = 1:nz

            I = squeeze(Imap(ix,iy,iz,:));

            % Convolve with total system response
            V = conv(I, h_total) * dt;

            % Truncate to nt
            Vmap(ix,iy,iz,:) = V(1:nt);

        end
    end
end

end