/time_derivs/d
/max_files/d
/out_params/d
/periodic/d
/courant=/d
/io_params/,/rubiks_params/{ /out_params/d }
/io_params/,/rubiks_params/{ :^/:d }
/atmosphere_/,/planet_/{ /mass_m_ea/d }
/planet_/,/patch_/{ /^ gamma/d }
