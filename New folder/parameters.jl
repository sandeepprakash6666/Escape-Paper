
module Bounds

    ##? Unscaled Initial POints/ Guesses
    x0_us = [30.0]
    z0_us = [60.0; 60.0; 60.0]
    u0_us = [0.0; 0.0]

    ##? Lower and Upper scales
    ls_x = [0.0]
    us_x = [100.0]

    ls_z = [0.0; 0.0; 0.0]
    us_z = [100.0; 100.0; 100.0]

    ls_u = [0.0; 0.0]
    us_u = [1.0; 10000000]

    #region-> Generic Code
    Nx = size(x0_us, 1)
    Nz = size(z0_us, 1)
    Nu = size(u0_us, 1)

    ##Scaled Points
    x0 = (x0_us - ls_x) ./ (us_x - ls_x)
    z0 = (z0_us - ls_z) ./ (us_z - ls_z)
    u0 = (u0_us - ls_u) ./ (us_u - ls_u)
    #endregion

    export  x0,      z0,     u0, 
            Nx,      Nz,     Nu, 
            ls_x,    ls_z,   ls_u, 
            us_x,    us_z,   us_u

end

module Model_parameters

    ## District Heat Parameters 
    ρ_dh = 1000
    Cp_dh = 4.18
    q_dh = 2

    ##Waste HEat Parameters
    ρ_wh = 1000
    Cp_wh = 4.18
    q_wh = 2

    T_dh_ret = 30
    T_dh_minSup = 60

    ##todo - make into arrays/ move out into varables
    Q_whb = vcat(1.2*ones(10,1), ones(10,1), 0.8*ones(10,1)) *1.2539989996092727e6
    V_tes = 100

    export  ρ_dh, Cp_dh, q_dh, 
            ρ_wh, Cp_wh, q_wh, 
            T_dh_ret, T_dh_minSup, Q_whb, V_tes

end