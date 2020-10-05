

function Var_bounds()

    ##? Unscaled Initial Points/ Guesses
    x0_us = [60.0;  60.0; 60.0]
    z0_us = [60.0]
    u0_us = [0.0;   0.0]

    des0_us = [100.0]

    ##? Lower and Upper scales
    ls_x = [0.0;    0.0;    0.0]
    us_x = [100.0;  100.0;  100.0]

    ls_z = [0.0]
    us_z = [100.0]

    ls_u = [0.0;    0.0]
    us_u = [1.0;    10000000]

    ls_des = [75.0]
    us_des = [125.0]

    #region-> #generic code (Collecting to Named Tuple)

        myTuple = (

        Nx   = size(x0_us, 1),
        Nz   = size(z0_us, 1),
        Nu   = size(u0_us, 1),
        Ndes = size(des0_us, 1),

        #scaling Bounds
        ls_x   = ls_x,
        us_x   = us_x, 
        ls_z   = ls_z, 
        us_z   = us_z, 
        ls_u   = ls_u,
        us_u   = us_u, 
        ls_des = ls_des, 
        us_des = us_des,

        ##Scaled Points
        x0   =  (x0_us - ls_x)      ./ (us_x - ls_x),
        z0   =  (z0_us - ls_z)      ./ (us_z - ls_z),
        u0   =  (u0_us - ls_u)      ./ (us_u - ls_u),
        des0 =  (des0_us - ls_des)  ./ (us_des - ls_des)

        )
    #endregion

    return myTuple

end

function Model_par()

    myTuple = (
        
    ## District Heat Parameters 
    ρ_dh = 1000,         
    Cp_dh = 4.18,        #kJ/(kg-K)
    q_dh = 10,            #m3/hr

    ##Waste Heat Parameters
    ρ_wh = 1000,
    Cp_wh = 4.18,
    q_wh = 10,

    T_dh_ret = 30,
    T_dh_minSup = 60,

    V_phb = 10.0,
    V_whb = 10.0
    )

    return myTuple

end

function Plant_par()

    myTuple = (
        
    ## District Heat Parameters 
    ρ_dh = 1000,         
    Cp_dh = 4.18,        #kJ/(kg-K)
    q_dh = 10,            #m3/hr

    ##Waste Heat Parameters
    ρ_wh = 1000,
    Cp_wh = 4.18,
    q_wh = 10,

    T_dh_ret = 30,
    T_dh_minSup = 60,

    V_phb = 10.0,
    V_whb = 10.0
    )

    return myTuple

end