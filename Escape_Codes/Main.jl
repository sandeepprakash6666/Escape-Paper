# using PlotlyJS
using Plots
include("OCP.jl")

dt = 1.0
NCP = 3
global N_Scen, N_Part = 1, 2


##* Create problem Objects
    #region-> #*Create Subproblem 1
    T01 = 0.0
    Tf1 = 30.0
    t_plot1 = collect(T01:dt:Tf1) 
    Q_whb1 = vcat(1.2*ones(10,1), 1.0*ones(10,1), 0.8*ones(10,1)) *1.2539999996092727e6
    
    P1 = Build_OCP(Q_whb1, Tf1 - T01 , (1,1))

    # Naming variables to access solution
    des1     = getindex(P1, :des)    
    x1       = getindex(P1, :x)
    x01      = getindex(P1, :x0)

    V_tes1   = getindex(P1, :V_tes)
    T_tes1   = getindex(P1, :T_tes)
    T_b1     = getindex(P1, :T_tes)
    T_phb1   = getindex(P1, :T_phb)
    T_whb1   = getindex(P1, :T_whb)
    α1       = getindex(P1, :α)
    Q_phb1   = getindex(P1, :Q_phb)
 
    #endregion

    #region-> #*Create Subproblem 2
    T02 = 30.0
    Tf2 = 60.0
    t_plot2 = collect(T02:dt:Tf2) 
    Q_whb2 = vcat(1.2*ones(10,1), 1.0*ones(10,1), 0.8*ones(10,1)) *1.2539999996092727e6 #todo - make into generic array or something similar

    P2 = Build_OCP(Q_whb2, Tf2 - T02, (1,2))

    # Naming variables to access solution
    des2    = getindex(P2, :des)    
    x2      = getindex(P2, :x)
    x02     = getindex(P2, :x0)

    V_tes2  = getindex(P2, :V_tes)
    T_tes2  = getindex(P2, :T_tes)
    T_b2    = getindex(P2, :T_tes)
    T_phb2  = getindex(P2, :T_phb)
    T_whb2  = getindex(P2, :T_whb)
    α2      = getindex(P2, :α)
    Q_phb2  = getindex(P2, :Q_phb)

    #endregion

    #region-> #* Create Central Problem
    T0 = T01
    Tf = Tf2

    Q_whb = vcat(Q_whb1, Q_whb2)

    t_plot = collect(T0:dt:Tf) 

    Centr = Build_OCP(Q_whb, Tf - T0, (1,1))
    # Naming variables to access solution
    des     = getindex(Centr, :des)    
    x       = getindex(Centr, :x)
    x0      = getindex(Centr, :x0)

    V_tes   = getindex(Centr, :V_tes)
    T_tes   = getindex(Centr, :T_tes)
    T_b     = getindex(Centr, :T_tes)
    T_phb   = getindex(Centr, :T_phb)
    T_whb   = getindex(Centr, :T_whb)
    α       = getindex(Centr, :α)
    Q_phb   = getindex(Centr, :Q_phb)
    #endregion

##* Solve Central Problem
    
    @NLobjective(Centr, Min, sum( Q_phb[nfe] for nfe in 1:60 ) + 40*(V_tes)^2)  #todo- generalize 60 and 40

    optimize!(Centr)
    JuMP.termination_status(Centr)
    JuMP.solve_time(Centr::Model)
    JuMP.objective_value(Centr)

    #extract solution to Julia variables
    star_x0 = JuMP.value.(x0) 
    star_x0_us = star_x0            .*  (100 - 0) .+ 0

    star_V_tes = JuMP.value(V_tes)
    star_T_tes = JuMP.value.(T_tes[:,NCP])
                star_T_tes = cat(star_x0_us[1], star_T_tes, dims = 1)    
    star_T_b    = JuMP.value.(T_b[:, NCP])
    star_T_phb  = JuMP.value.(T_phb[:, NCP])
    star_T_whb  = JuMP.value.(T_whb[:, NCP])

    star_α      = JuMP.value.(α[:])
    star_Q_phb  = JuMP.value.(Q_phb[:])


    #plot Central solution
    plotly()
    plot(t_plot, star_T_tes, ylim = [55,70])

##* Solve subproblems

z_des_us = [120.0]
z_diff_us = [60.0]

z_des  = (z_des_us - ls_des) ./ (us_des - ls_des)
z_diff = (z_diff_us - ls_x)  ./ (us_x   - ls_x  )
lambda1_des, lambda2_des    = 0.0.*z_des,  0.0.*z_des
lambda1_diff, lambda2_diff  = 0.0.*z_diff, 0.0.*z_diff

rho = 1e6

    #region-> declaring arrays for Plotting
    plot_z_des = [z_des]
    plot_z_diff = [z_diff]
    plot_lambda1_des, plot_lambda2_des = [lambda1_des], [lambda2_des]
    plot_lambda1_diff, plot_lambda2_diff = [lambda1_diff], [lambda2_diff]

    plot_V_tes1 = [120.0]
    plot_V_tes2 = [120.0]
    plot_T_tes1 = NaN*t_plot1
    plot_T_tes2 = NaN*t_plot2
    #endregion

##* ADMM Iterations
NIter = 10
for ADMM_k = 1:NIter
    global z_des, z_diff
    global lambda1_des, lambda2_des
    global lambda1_diff,lambda2_diff
    global plot_V_tes1, plot_V_tes2
    global plot_T_tes1, plot_T_tes2
    
    #region-> #*Solve SP 1
    
    @NLobjective(P1, Min, sum( Q_phb1[nfe] for nfe in 1:30 ) + 20*(V_tes1)^2    + lambda1_des[1]*   (des1[1] - z_des[1])                + rho/2*(des1[1] - z_des[1])^2  
                                                                                + lambda1_diff[1]*  (x1[1, end, end] - z_diff[1])       + rho/2*(x1[1, end, end] - z_diff[1])^2  )     #todo - Make the penalty variables between 0-1
    
    optimize!(P1)
    JuMP.termination_status(P1)
    JuMP.solve_time(P1::Model)
    JuMP.objective_value(P1)
        
    ## extract solution to Julia variables
        #scaled values
        star_des1 = JuMP.value.(des1)
        star_x01 = JuMP.value.(x01) 
        star_x1 = JuMP.value.(x1)
        
        #unscaled values
        star_x01_us = star_x01            .*  (100 - 0) .+ 0
        star_V_tes1 = JuMP.value(V_tes1)
        star_T_tes1 = JuMP.value.(T_tes1[:,NCP])
                    star_T_tes1 = cat(star_x01_us[1], star_T_tes1, dims = 1)     
        star_T_b1    = JuMP.value.(T_b1[:, NCP])
        star_T_phb1  = JuMP.value.(T_phb1[:, NCP])
        star_T_whb1  = JuMP.value.(T_whb1[:, NCP])

        star_α1      = JuMP.value.(α1[:])
        star_Q_phb1  = JuMP.value.(Q_phb1[:])
    #endregion

    #region-> #*Solve SP2 

    @NLobjective(P2, Min, sum( Q_phb2[nfe] for nfe in 1:30 ) + 20*(V_tes2)^2    + lambda2_des[1] *  (des2[1] - z_des[1])                + rho/2*(des2[1] - z_des[1])^2  
                                                                                + lambda2_diff[1]*  (x02[1]  - z_diff[1])               + rho/2*(x02[1] - z_diff[1])^2  )   #?1st element linked here i.e. x0

    optimize!(P2)
    JuMP.termination_status(P2)
    JuMP.solve_time(P2::Model)
    JuMP.objective_value(P2)

    ## extract solution to Julia variables
        #scaled values
        star_des2 = JuMP.value.(des2)
        star_x02 = JuMP.value.(x02) 
        star_x2 = JuMP.value.(x2)
        
        #unscaled Values
        star_x02_us = star_x02            .*  (100.0 - 0.0) .+ 0.0          #todo - scaling to be made automatic
        star_V_tes2 = JuMP.value(V_tes2)
        star_T_tes2 = JuMP.value.(T_tes2[:,NCP])
                    star_T_tes2 = cat(star_x02_us[1], star_T_tes2, dims = 1)     
        star_T_b2    = JuMP.value.(T_b2[:, NCP])
        star_T_phb2  = JuMP.value.(T_phb2[:, NCP])
        star_T_whb2  = JuMP.value.(T_whb2[:, NCP])

        star_α2      = JuMP.value.(α2[:])
        star_Q_phb2  = JuMP.value.(Q_phb2[:])
    #endregion

    #* ADMM Updates
    z_des[1]    = (star_des1[1] + star_des2[1])/2
    z_diff[1]   = (star_x1[1, end, end] + star_x02[1])/2

    lambda1_des[1] = lambda1_des[1] + rho*(star_des1[1] - z_des[1])
    lambda2_des[1] = lambda2_des[1] + rho*(star_des2[1] - z_des[1])

    lambda1_diff[1] = lambda1_diff[1] + rho*(star_x1[1, end, end]   - z_diff[1])
    lambda2_diff[1] = lambda2_diff[1] + rho*(star_x02[1]            - z_diff[1])

        #region-> Storing in Plots
        plot_V_tes1 = append!(plot_V_tes1, star_V_tes1)
        plot_V_tes2 = append!(plot_V_tes2, star_V_tes2)
        
        plot_T_tes1 = cat(plot_T_tes1, star_T_tes1, dims = 2)
        plot_T_tes2 = cat(plot_T_tes2, star_T_tes2, dims = 2)
        # append!(plot_lambda1, lambda1)
        # append!(plot_lambda2, lambda2)
        # append!(plot_z, z)
        #endregion

        #todo - Setting warm start - How to?
end


##* Plotting ADMM iterations - Summary
    #choose backend for plots
    plotly()
    # gr()

    #final point
    p11 = plot(t_plot1,  plot_T_tes1[:,end], ylim = [55,70] )
    p11 = plot!(t_plot2, plot_T_tes2[:,end], ylim = [55,70])

    plot_V_tes1
    plot_V_tes2
    plot_T_tes1
    plot_T_tes2

    plot_Iter = collect(1:NIter)
    p21 = plot(plot_Iter, plot_V_tes1[:])
    # gui(fig1)

    # plot_lambda1
    # plot_z

##* Plotting Complete Profiles

    #region #*Plot Subproblem 1
    # #Differential States
    # p21 = plot(t_plot1, star_T_tes1)
    # #Manipulated Variables
    # p22 = plot(t_plot1[1:end-1], star_α1,
    #                                                                     marker = true, linetype = :steppost)
    # p23 = plot(t_plot1[1:end-1], star_Q_phb1,        
    #                                                                     marker = true, linetype = :steppost)
    # #Algebraic States
    # p24 = plot(t_plot1[2:end], star_T_b1, 
    #                                                                     marker = true)
    # p24 = plot!(t_plot1[2:end], star_T_phb1) 
    #                                                                     #marker = false)
    # p25 = plot(t_plot1[2:end], star_T_whb1) 
    #                                                             # marker = false)

    # fig_P2 = plot(p21, p22, p23, p24, p25, layout = (5, 1))
    # # gui(fig1)
    #endregion

    #region #*Plot Subproblem 2
    #     #Differential States
    #     p31 = plot(t_plot2, star_T_tes2                  )
    #     #Manipulated Variables
    #     p32 = plot(t_plot2[1:end-1], star_α2,
    #                                                                         marker = true, linetype = :steppost)
    #     p33 = plot(t_plot2[1:end-1], star_Q_phb2,
    #                                                                         marker = true, linetype = :steppost)
    #     #Algebraic States
    #     p34 = plot(t_plot2[2:end], star_T_b2, 
    #                                                                         marker = true)
    #     p34 = plot!(t_plot2[2:end], star_T_phb2) 
    #                                                                         #marker = false)
    #     p35 = plot(t_plot2[2:end], star_T_whb2) 
    #                                                                 # marker = false)

    # fig_P3 = plot(p31, p32, p33, p34, p35, layout = (5, 1))


    # fig_3 = plot(fig_P1, fig_P2, layout = (1,2))

    # fig_4 = plot(p11, p21)
    #endregion

