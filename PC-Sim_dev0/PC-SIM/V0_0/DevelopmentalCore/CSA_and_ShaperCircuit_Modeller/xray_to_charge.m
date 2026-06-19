function Q_event = xray_to_charge(E_keV, w_eV)
e     = 1.602e-19;

N_e = (E_keV*1e3) / w_eV;
Q_event = N_e * e;

end
