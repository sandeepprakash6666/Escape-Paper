# using PlotlyJS
using Plots
include("OCP.jl")

dt = 1.0
NCP = 3
global N_Scen, N_Part = 1, 2
Obj_scaling = 2e1

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
    u1       = getindex(P1, :u)

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
    Q_whb2 = vcat(1.3*ones(10,1), 1.0*ones(10,1), 0.7*ones(10,1)) *1.2539999996092727e6 #todo - make into generic array or something similar

    P2 = Build_OCP(Q_whb2, Tf2 - T02, (1,2))

    # Naming variables to access solution
    des2    = getindex(P2, :des)    
    x2      = getindex(P2, :x)
    x02     = getindex(P2, :x0)
    u2      = getindex(P2, :u)

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
    u       = getindex(Centr, :u) 
    
    V_tes   = getindex(Centr, :V_tes)
    T_tes   = getindex(Centr, :T_tes)
    T_b     = getindex(Centr, :T_tes)
    T_phb   = getindex(Centr, :T_phb)
    T_whb   = getindex(Centr, :T_whb)
    α       = getindex(Centr, :α)
    Q_phb   = getindex(Centr, :Q_phb)
    #endregion

##* Solve Central Problem
    
    @NLobjective(Centr, Min, Obj_scaling*(sum( u[2,nfe] for nfe in 1:60 ) + 0.1*(des[1])^2) )  #todo- generalize

    optimize!(Centr)
    JuMP.termination_status(Centr)
    JuMP.solve_time(Centr::Model)
    star_Obj = JuMP.objective_value(Centr)

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
    plot(t_plot, star_T_tes, ylim = [55,70], label = "Ttes - Centr", title = "Diff state, Central" )

##* Solve subproblems using ADMM
rho = 2.0*Obj_scaling

z_des_us = [100.0]
z_diff_us = [58.33155]

    z_des  = (z_des_us - ls_des) ./ (us_des - ls_des)
    z_diff = (z_diff_us - ls_x)  ./ (us_x   - ls_x  )

μ1_des   = 0.0.*z_des
μ1_diff  = 0.0.*z_diff

    μ2_des  = 0.0.*z_des
    μ2_diff = 0.0.*z_diff

    #region-> declaring arrays for Plotting
    plot_rho = copy([rho])

    plot_z_des  = copy(z_des)
    plot_z_diff = copy(z_diff)

    plot_des1 = copy(z_des)
    plot_x1 = NaN*t_plot1[2:end]
    plot_x01 = copy(z_diff)
    plot_μ1_des = copy(μ1_des)
    plot_μ1_diff = copy(μ1_diff)
    plot_Obj1 = [NaN]
        
        plot_des2 = copy(z_des)
        plot_x2 = NaN*t_plot2[2:end]
        plot_x02 = copy(z_diff)
        plot_μ2_des = copy(μ2_des)
        plot_μ2_diff = copy(μ2_diff)
        plot_Obj2 = [NaN]

    plot_V_tes1 = copy(z_des_us)
    plot_T_tes1 = NaN*t_plot1

        plot_V_tes2 = copy(z_des_us)
        plot_T_tes2 = NaN*t_plot2
    #endregion

##* ADMM Iterations
NIter = 50
Tot_time_in_ADMM = @elapsed for ADMM_k = 2:NIter
    #ADMM values
    global rho
    global plot_rho
    global z_des, z_diff
    global μ1_des, μ2_des
    global μ1_diff,μ2_diff
        #Plotting - scaled
        global plot_z_des, plot_z_diff
        global plot_des1, plot_des2
        global plot_x01, plot_x02
        global plot_x1, plot_x2
        global plot_μ1_des, plot_μ2_des
        global plot_μ1_diff, plot_μ2_diff
        global plot_Obj1, plot_Obj2
        #Plotting - Unscaled
        global plot_V_tes1, plot_V_tes2
        global plot_T_tes1, plot_T_tes2
    
    #region-> #*Solve SP 1
    
    @NLobjective(P1, Min, Obj_scaling*(sum( u1[2,nfe] for nfe in 1:30 ) + 0.05*(des1[1])^2)     + μ1_des[1]*   (des1[1] - z_des[1])                + rho/2*(des1[1] - z_des[1])^2  
                                                                                                + μ1_diff[1]*  (x1[1, end, end] - z_diff[1])       + rho/2*(x1[1, end, end] - z_diff[1])^2  )     #todo- generalize
    
    optimize!(P1)
    JuMP.termination_status(P1)
    JuMP.solve_time(P1::Model)
    star_Obj1 = JuMP.objective_value(P1)
        
    ## extract solution to Julia variables
        #scaled values
        star_des1 = JuMP.value.(des1)
        star_x01 = JuMP.value.(x01) 
        star_x1 = JuMP.value.(x1)[:,:,NCP]
        
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

    @NLobjective(P2, Min, Obj_scaling*(sum( u2[2,nfe] for nfe in 1:30 ) + 0.05*(des2[1])^2)     + μ2_des[1] *  (des2[1] - z_des[1])                + rho/2*(des2[1] - z_des[1])^2  
                                                                                                + μ2_diff[1]*  (x02[1]  - z_diff[1])               + rho/2*(x02[1] - z_diff[1])^2  )   #?1st element linked here i.e. x0

    optimize!(P2)
    JuMP.termination_status(P2)
    JuMP.solve_time(P2::Model)
    star_Obj2 = JuMP.objective_value(P2)

    ## extract solution to Julia variables
        #scaled values
        star_des2 = JuMP.value.(des2)
        star_x02 = JuMP.value.(x02) 
        star_x2 = JuMP.value.(x2)[:,:,NCP]
        
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


    μ1_des[1] = μ1_des[1] + rho*(star_des1[1] - z_des[1])
    μ2_des[1] = μ2_des[1] + rho*(star_des2[1] - z_des[1])

    μ1_diff[1] = μ1_diff[1] + rho*(star_x1[1, end, end]   - z_diff[1])
    μ2_diff[1] = μ2_diff[1] + rho*(star_x02[1]            - z_diff[1])

    prim_res_des  = star_des1[1]  - z_des[1] 
    prim_res_diff = star_x1[1, end, end] - z_diff[1]
    prim_res =  (prim_res_des[1]^2 + prim_res_diff[1]^2)^0.5

    dual_res_des  = rho*(z_des[1]  - plot_z_des[end])
    dual_res_diff = rho*(z_diff[1] - plot_z_diff[end])
    dual_res = (dual_res_des[1]^2 + dual_res_diff[1]^2)^0.5

    ##* rho update heuristic
    if prim_res > 10*dual_res
        rho = rho*2
    elseif dual_res > 10*prim_res
        rho = rho/2
    else
    end
    append!(plot_rho, rho)

        #region-> #*Storing in Plots
            
            #scaled values
            append!(plot_z_des, z_des)
            append!(plot_z_diff, z_diff)
            
            append!(plot_des1, star_des1)
            append!(plot_x01, star_x01)
            plot_x1 = cat(plot_x1, transpose(star_x1), dims = 2)
            append!(plot_μ1_des, μ1_des)
            append!(plot_μ1_diff, μ1_diff)
            append!(plot_Obj1, star_Obj1)

                append!(plot_des2, star_des2)
                append!(plot_x02, star_x02)
                plot_x2 = cat(plot_x2, transpose(star_x2), dims = 2)
                append!(plot_μ2_des, μ2_des)
                append!(plot_μ2_diff, μ2_diff)
                append!(plot_Obj2, star_Obj2)

            #Unscaled values
            append!(plot_V_tes1, star_V_tes1)
            plot_T_tes1 = cat(plot_T_tes1, star_T_tes1, dims = 2)

                append!(plot_V_tes2, star_V_tes2)
                plot_T_tes2 = cat(plot_T_tes2, star_T_tes2, dims = 2)

        #endregion

        #region-> #*Warm Starting

            # for ndes in 1:Ndes, nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:30, ncp in 1:3 #todo- generalise NFE and NCP
            # set_start_value(des1[ndes],                star_des1[ndes])
            # set_start_value(x01[nx],                   star_x01[nx])
            # set_start_value(x1[nx, nfe, ncp],          star_x1[nx, nfe])
            # # set_start_value(z1[nz, nfe, ncp],          z_guess[nz])
            # # set_start_value(u1[nu, nfe],               u_guess[nu])

            #         set_start_value(des2[ndes],                star_des2[ndes])
            #         set_start_value(x02[nx],                   star_x02[nx])
            #         set_start_value(x2[nx, nfe, ncp],          star_x2[nx, nfe])
            #         # set_start_value(z2[nz, nfe, ncp],          z_guess[nz])
            #         # set_start_value(u2[nu, nfe],               u_guess[nu])
            # end

        #endregion

end
plot_rho
##* Calculating ADMM - Summary Stats #todo - add here
    
                                #region->
                                plot_μ1_des
                                plot_μ1_diff
                                #endregion
    #Penalty
        plot_Aug_terms1   = NaN*zeros(NIter)
        plot_Aug_terms2   = NaN*zeros(NIter)
        for nIter in 1:NIter   
            # nIter = 2
            plot_Aug_terms1[nIter] =  plot_μ1_des[nIter] *(plot_des1[nIter]   - plot_z_des[nIter])       + rho/2*(plot_des1[nIter]   - plot_z_des[nIter])^2   + 
                                      plot_μ1_diff[nIter]*(plot_x1[end,nIter] - plot_z_diff[nIter])      + rho/2*(plot_x1[end,nIter] - plot_z_diff[nIter])^2

            plot_Aug_terms2[nIter] =  plot_μ2_des[nIter] *(plot_des2[nIter] - plot_z_des[nIter])         + rho/2*(plot_des2[nIter] - plot_z_des[nIter])^2     + 
                                      plot_μ2_diff[nIter]*(plot_x02[nIter]  - plot_z_diff[nIter])        + rho/2*(plot_x02[nIter]  - plot_z_diff[nIter])^2
        end 

    #Original Objective function value
        plot_f1_star    = NaN*zeros(NIter)
        plot_f2_star    = NaN*zeros(NIter)
        for nIter in 1:NIter  
            # nIter = 2 
            plot_f1_star[nIter] = plot_Obj1[nIter] - plot_Aug_terms1[nIter] 
            plot_f2_star[nIter] = plot_Obj2[nIter] - plot_Aug_terms2[nIter] 
        end
    
    #Primal Residual
        plot_prim_res_des = NaN*zeros(NIter)
        plot_prim_res_diff = NaN*zeros(NIter)
        plot_prim_res = NaN*zeros(NIter)
        for nIter in 2:NIter
            # nIter = 2
            plot_prim_res_des[nIter]  = plot_des1[nIter]   - plot_z_des[nIter]
            plot_prim_res_diff[nIter] = plot_x1[end,nIter] - plot_z_diff[nIter]

            plot_prim_res[nIter] = (plot_prim_res_des[nIter]^2 + plot_prim_res_diff[nIter]^2)^0.5
        end
    
    #Dual Residual
        plot_dual_res_des = NaN*zeros(NIter)
        plot_dual_res_diff = NaN*zeros(NIter)
        plot_dual_res = NaN*zeros(NIter)
        for nIter in 2:NIter
            # nIter = 2
            plot_dual_res_des[nIter]  = plot_rho[nIter]*(plot_z_des[nIter]  - plot_z_des[nIter-1])
            plot_dual_res_diff[nIter] = plot_rho[nIter]*(plot_z_diff[nIter] - plot_z_diff[nIter-1])
            
            plot_dual_res[nIter] = (plot_dual_res_des[nIter]^2 + plot_dual_res_diff[nIter]^2)^0.5
        end
                    
                    plot_Aug_terms1
                    plot_Aug_terms1
                    
                    plot_f1_star
                    plot_f2_star

                    plot_prim_res_des
                    plot_prim_res_diff
                    plot_prim_res
                    plot_dual_res_des
                    plot_dual_res_diff
                    plot_dual_res

##* Plot ADMM itreations Summary

    #choose backend for plots
    plotly()
    # gr()

    #Profiles - Differential States
        p11 = plot(t_plot1,  plot_T_tes1[:,end], ylim = [55,70], label = "Ttes - SP1" )
        p11 = plot!(t_plot2, plot_T_tes2[:,end], ylim = [55,70], label = "Ttes - SP2", title = "Diff state, rho = $rho" )

    #Iteration Profiles - Design Variable
        plot_Iter = collect(0:NIter-1)
        p21 = plot(plot_Iter, plot_V_tes1[:],               label = "des - SP1")
        p21 = plot!(plot_Iter, plot_V_tes2[:],              label = "des SP1")
        p21 = plot!(plot_Iter, star_V_tes.*ones(NIter),     label = "des Centr", title = "Des var, rho = $rho" )

    #Primal and Dual Infeasibility
        p31 = plot(plot_Iter, plot_prim_res,                label  = "primal infeasibility", yscale = :log10)
        p31 = plot!(plot_Iter, plot_dual_res,               label = "dual infeasibility",    yscale = :log10, title = "Infeasibility, rho = $rho")

    #Objective Function Values
        p41 = plot(plot_Iter, plot_f1_star .+ plot_f2_star, label = "f - Distribuited")
        p41 = plot!(plot_Iter, star_Obj.*ones(NIter),       label = "f - Central Sol'n" )
        p41 = plot!(plot_Iter, plot_Obj1 .+ plot_Obj2,      label = "Augmented Lagrange", title = "Objective, rho = $rho")
    
        p51 = plot(plot_Iter, plot_f1_star,                 label ="f1 star")
        p51 = plot!(plot_Iter, plot_f2_star,                label ="f2 star")
        p51 = plot!(plot_Iter, plot_Aug_terms1,             label = "f1 Aug Terms")
        p51 = plot!(plot_Iter, plot_Aug_terms2,             label = "f2 Aug Terms", title = "SP1 objective, rho = $rho")

        p61 = plot(plot_Iter, plot_μ1_des,                  label = "μ1_des")
        p61 = plot!(plot_Iter, plot_μ2_des,                 label = "μ2_des")
        p61 = plot!(plot_Iter, plot_μ1_diff,                label = "μ1 diff")
        p61 = plot!(plot_Iter, plot_μ2_diff,                label = "μ2 diff", title = "multipliers, rho = $rho")

        p71 = plot(plot_Iter, plot_rho,                     title = "Rho updates")
##* Display Plots
p11
p21
p31
p41
p51
p61
display(Tot_time_in_ADMM)


