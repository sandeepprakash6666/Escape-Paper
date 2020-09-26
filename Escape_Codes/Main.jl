# using PlotlyJS
dt = 1.0
NCP = 3
using Plots
include("OCP.jl")


global NS, NP = 1, 2
Obj_scaling = 1e0

##* Create Model Objects
    Tf = 60
    NFE = convert(Int32, (Tf - 0)/dt)

    #make Profiles
    plot_t = collect(0: dt: Tf)
    Q_whb = hcat(1.2*ones(1,10), 1.0*ones(1,10), 0.8*ones(1,10),    1.3*ones(1,10), 1.0*ones(1,10), 0.7*ones(1,10)) *1.2539999996092727e6 

    #Make central Problem 
    Centr = Build_OCP(Q_whb, plot_t[end] - plot_t[1] , (1,1))

        #region -> #*make Dictionary of Subproblem Models
            SP = Dict()
            SP_len = convert(Int32,size(Q_whb)[2]/ NP)
            for nS in 1:NS, nP in 1:NP
            # nS = 1
            # nP = 3
                Q_whb_nS_nP = Q_whb[nS, SP_len*(nP-1)+1 : SP_len*(nP) ] 
                t_nS_nP = plot_t[SP_len*(nP-1)+1 : SP_len*(nP)+1 ] 

                SP[(nS,nP)] = Build_OCP(Q_whb_nS_nP, t_nS_nP[end] - t_nS_nP[1] , (nS,nP))
            end
                            SP
                            SP[(1,1)]
                            SP[(1,2)]
                            # SP[(1,3)]
        #endregion

        #region-> #*Create Named references

        #Subproblems
            SP_des  = Dict()
            SP_x    = Dict()
            SP_x0   = Dict()
            SP_u    = Dict()

            SP_V_tes = Dict()
            SP_T_tes = Dict()
            for nS in 1:NS, nP in 1:NP
            # nS = 1
            # nP = 1
                #Scaled variables
                SP_des[(nS,nP)] = getindex(SP[(nS,nP)], :des) 
                SP_x[(nS,nP)]   = getindex(SP[(nS,nP)], :x)
                SP_x0[(nS,nP)]  = getindex(SP[(nS,nP)], :x0)
                SP_u[(nS,nP)]   = getindex(SP[(nS,nP)], :u)

                #unscaled variables
                SP_V_tes[(nS,nP)]   = getindex(SP[(nS,nP)], :V_tes)
                SP_T_tes[(nS,nP)]   = getindex(SP[(nS,nP)], :T_tes)[:,end]
            end
                            
                            SP_des[(1,1)][1]
                            SP_x[(1,1)][1,end,end] 

        #Central Solver
                #Scaled variables
                Centr_des = getindex(Centr, :des) 
                Centr_x   = getindex(Centr, :x)
                Centr_x0  = getindex(Centr, :x0)
                Centr_u   = getindex(Centr, :u)
        
                #unscaled variables
                Centr_V_tes   = getindex(Centr, :V_tes)
                Centr_T_tes   = getindex(Centr, :T_tes)[:,end]

        #endregion

##* Solve Central Problem)

    @NLobjective(Centr, Min, Obj_scaling*(sum( Centr_u[2,nfe] for nfe in 1:Tf ) + 0.1*(Centr_des[1])^2) )

    optimize!(Centr)
    JuMP.termination_status(Centr)
    JuMP.solve_time(Centr::Model)
    star_Obj = JuMP.objective_value(Centr)

    #extract solution to Julia variables
    star_x0 = JuMP.value.(Centr_x0) 
    star_x0_us = star_x0            .*  (100 - 0) .+ 0

    star_V_tes = JuMP.value(Centr_V_tes)
    star_T_tes = JuMP.value.(Centr_T_tes)
                 star_T_tes = cat(star_x0_us[1], star_T_tes, dims = 1)    
    
    #plot Central solution
    println("Central Solution is $star_V_tes")

    plotly()
    plot(plot_t, star_T_tes, ylim = [55,70], label = "Ttes - Centr", title = "Diff state, Central" )

##* Solve subproblems using ADMM
rho = 2.0*Obj_scaling

z_des_us = [100.0]
z_diff_us = [60.0]

    z_des = NaN.*zeros(Ndes)
    z_diff = NaN.*zeros(NS, NP-1, Nx)
    for ndes in 1:Ndes, nx in 1:Nx
        z_des[ndes]         = (z_des_us[ndes] - ls_des[ndes]) / (us_des[ndes] - ls_des[ndes])
        z_diff[:, :, nx]   .= (z_diff_us[nx]  - ls_x[nx])     / (us_x[nx]     - ls_x[nx]  )
    end

    μ_des    = zeros(NS, NP, Ndes)
    μ_diff_L = zeros(NS, NP, Nx)
    μ_diff_R = zeros(NS, NP, Nx)

    #region-> Initializing arrays for Plotting
        plot_rho = []
        plot_z_des  = []
        plot_z_diff = []

        #scaled variables
        plot_Obj = []
        plot_des = []
        plot_x = []
        plot_x0 = []
        plot_u = []

        #Unscaled Variables
        plot_V_tes = []
        plot_T_tes = []

        #Multipliers
        plot_μ_des = []
        plot_μ_diff_L = []
        plot_μ_diff_R = []
        
        push!(plot_rho, rho)
        push!(plot_z_des,  copy(z_des) )
        push!(plot_z_diff, copy(z_diff))

    #endregion

##* ADMM Iterations

NIter = 100
Tot_time_in_ADMM = @elapsed for ADMM_k = 1:NIter
    #ADMM values
    global rho
    global z_des, z_diff
    global μ_des, μ_diff_L, μ_diff_R
    
    #region-> #*Solve all partition problems
        star_SP_Obj = zeros(NS, NP)*NaN
        star_SP_des = zeros(NS, NP, Ndes)*NaN
        star_SP_x0  = zeros(NS, NP, Nx)*NaN
        star_SP_x   = zeros(NS, NP, Nx, convert(Int32,NFE/NP))*NaN
        star_SP_u   = zeros(NS, NP, Nu, convert(Int32,NFE/NP))*NaN

        star_SP_V_tes = zeros(NS, NP)
        star_SP_T_tes = zeros(NS, NP, convert(Int32,NFE/NP))*NaN

        for nS in 1:NS, nP in 1:NP
                
            if nP == 1          #* First Partition
                # nS = 1 
                # nP = 1
                @NLobjective(SP[(nS,nP)], Min,  
                                
                                Obj_scaling*(sum(SP_u[(nS,nP)][2,nfe]  for nfe in 1:SP_len)                   + (0.1/(NS*NP))*(SP_des[(nS,nP)][1])^2)
                                                
                                + sum(  μ_des[nS, nP, ndes] *(SP_des[(nS,nP)][ndes]     - z_des[ndes])        + rho*(SP_des[(nS,nP)][ndes]      - z_des[ndes])^2             for ndes in 1:Ndes)
                                
                                + sum(  μ_diff_R[nS, nP, nx]*(SP_x[(nS,nP)][nx,end,end] - z_diff[nS,nP,nx])   + rho*(SP_x[(nS,nP)][nx,end,end]  - z_diff[nS,nP,nx])^2        for nx in 1:Nx) 
                            )
                
            elseif nP == NP     #* Last Partition
                # nS = 1 
                # nP = 3
                @NLobjective(SP[(nS,nP)], Min,  
                                
                                Obj_scaling*(sum(SP_u[(nS,nP)][2,nfe]  for nfe in 1:SP_len)                   + (0.1/(NS*NP))*(SP_des[(nS,nP)][1])^2)
                                                
                                + sum(  μ_des[nS, nP, ndes] *(SP_des[(nS,nP)][ndes] - z_des[ndes])          + rho*(SP_des[(nS,nP)][ndes]    - z_des[ndes])^2            for ndes in 1:Ndes)
                                
                                + sum(  μ_diff_L[nS, nP, nx]*(SP_x0[(nS,nP)][nx]    - z_diff[nS,nP-1,nx])   + rho*(SP_x0[(nS,nP)][nx]       - z_diff[nS,nP-1,nx])^2     for nx in 1:Nx) 
                            )

            else                #* Between Partitions
                # nS = 1 
                # nP = 2
                @NLobjective(SP[(nS,nP)], Min,  
                                
                                Obj_scaling*(sum(SP_u[(nS,nP)][2,nfe]  for nfe in 1:SP_len)                   + (0.1/(NS*NP))*(SP_des[(nS,nP)][1])^2)
                                                
                                + sum(  μ_des[nS, nP, ndes] *(SP_des[(nS,nP)][ndes]     - z_des[ndes])          + rho*(SP_des[(nS,nP)][ndes]      - z_des[ndes])^2            for ndes in 1:Ndes)
                                
                                + sum(  μ_diff_L[nS, nP, nx]*(SP_x0[(nS,nP)][nx]        - z_diff[nS,nP-1,nx])   + rho*(SP_x0[(nS,nP)][nx]         - z_diff[nS,nP-1,nx])^2     for nx in 1:Nx) 
                                
                                + sum(  μ_diff_R[nS, nP, nx]*(SP_x[(nS,nP)][nx,end,end] - z_diff[nS,nP,nx])     + rho*(SP_x[(nS,nP)][nx,end,end]  - z_diff[nS,nP,nx])^2       for nx in 1:Nx) 
                            )
                
            end
            
            #solve problem
            optimize!(SP[(nS,nP)])
            JuMP.termination_status(SP[(nS,nP)])
            JuMP.solve_time(SP[(nS,nP)]::Model)
            star_SP_Obj[nS,nP] = JuMP.objective_value(SP[(nS,nP)])

            #solution from SP
            star_SP_des[nS, nP, :]  = JuMP.value.(SP_des[(nS,nP)])
            star_SP_x0[nS, nP, :]   = JuMP.value.(SP_x0[(nS,nP)])
            star_SP_x[nS, nP, :, :] = JuMP.value.(SP_x[(nS,nP)])[:,:,NCP]
            star_SP_u[nS, nP, :, :] = JuMP.value.(SP_u[(nS,nP)])
        
            star_SP_V_tes[nS,nP]    = JuMP.value(SP_V_tes[(nS,nP)])
            star_SP_T_tes[nS,nP,:]  = JuMP.value.(SP_T_tes[(nS,nP)])
            
        end
    #endregion

        star_SP_Obj
        star_SP_des
        star_SP_x0
        star_SP_x
        star_SP_u

        star_SP_V_tes
        star_SP_T_tes

    #region->Oldcode

        # @NLobjective(P1, Min, Obj_scaling*(sum( u1[2,nfe] for nfe in 1:20 ) + 0.1/3*(des1[1])^2)     + μ1_des[1]*   (des1[1] - z_des[1])            + rho/2*(des1[1] - z_des[1])^2  
        #                                                                                              + μ1_diff[1]*  (x1[1, end, end] - z_diff[1])   + rho/2*(x1[1, end, end] - z_diff[1])^2  )     #todo- generalize
        
        # optimize!(P1)
        # JuMP.termination_status(P1)
        # JuMP.solve_time(P1::Model)
        # star_Obj1 = JuMP.objective_value(P1)
            
        # ## extract solution to Julia variables
        #     #scaled values
        #     star_des1 = JuMP.value.(des1)
        #     star_x01 = JuMP.value.(x01) 
        #     star_x1 = JuMP.value.(x1)[:,:,NCP]
            
        #     #unscaled values
        #     star_x01_us = star_x01            .*  (100 - 0) .+ 0
        #     star_V_tes1 = JuMP.value(V_tes1)
        #     star_T_tes1 = JuMP.value.(T_tes1[:,NCP])
        #                 star_T_tes1 = cat(star_x01_us[1], star_T_tes1, dims = 1)     
        #     star_T_b1    = JuMP.value.(T_b1[:, NCP])
        #     star_T_phb1  = JuMP.value.(T_phb1[:, NCP])
        #     star_T_whb1  = JuMP.value.(T_whb1[:, NCP])

        #     star_α1      = JuMP.value.(α1[:])
        #     star_Q_phb1  = JuMP.value.(Q_phb1[:])

        # ##Solve SP2 

        # @NLobjective(P2, Min, Obj_scaling*(sum( u2[2,nfe] for nfe in 1:20 ) + 0.1/3*(des2[1])^2)    + μ2_des[1] *  (des2[1] - z_des[1])             + rho/2*(des2[1] - z_des[1])^2  
                                                                                                    
        #                                                                                             + μ2_diff[1]*  (x02[1]  - z_diff[1])            + rho/2*(x02[1] - z_diff[1])^2  
        #                                                                                             + μ2_diff[2]*  (x2[1, end, end]  - z_diff[2])   + rho/2*(x2[1, end, end]  - z_diff[2])^2 )   #?1st element linked here i.e. x0

        # optimize!(P2)
        # JuMP.termination_status(P2)
        # JuMP.solve_time(P2::Model)
        # star_Obj2 = JuMP.objective_value(P2)

        # ## extract solution to Julia variables
        #     #scaled values
        #     star_des2 = JuMP.value.(des2)
        #     star_x02 = JuMP.value.(x02) 
        #     star_x2 = JuMP.value.(x2)[:,:,NCP]
            
        #     #unscaled Values
        #     star_x02_us = star_x02            .*  (100.0 - 0.0) .+ 0.0         
        #     star_V_tes2 = JuMP.value(V_tes2)
        #     star_T_tes2 = JuMP.value.(T_tes2[:,NCP])
        #                 star_T_tes2 = cat(star_x02_us[1], star_T_tes2, dims = 1)     
        #     star_T_b2    = JuMP.value.(T_b2[:, NCP])
        #     star_T_phb2  = JuMP.value.(T_phb2[:, NCP])
        #     star_T_whb2  = JuMP.value.(T_whb2[:, NCP])

        #     star_α2      = JuMP.value.(α2[:])
        #     star_Q_phb2  = JuMP.value.(Q_phb2[:])
    #endregion
 
    #* ADMM Updates
    #region ->#*z-Updates  
        #Design
        for ndes in 1:Ndes
            # ndes = 1
            z_des[ndes] = sum(star_SP_des[nS, nP, ndes] for nS in 1:NS, nP in 1:NP )/(NS*NP)
        end

        #Differential
        for nS in 1:NS, nP in 1:NP-1, nx in 1:Nx
        # nS = 1; nP = 1; nx = 1
            z_diff[nS,nP,nx] = (star_SP_x[nS, nP, nx, end]  +   star_SP_x0[nS, nP+1, nx] )/2   
        end
        
    z_des
    z_diff
    #endregion

    #region ->#*Primal Residuals
        #Design
        prim_res_des    = zeros(NS, NP, Ndes)
        for nS in 1:NS, nP in 1:NP, ndes in 1:Ndes
            prim_res_des[nS, nP, ndes] = star_SP_des[nS, nP, ndes] - z_des[ndes]
        end
        
        #Differential
        prim_res_diff_L = zeros(NS, NP, Nx)
        prim_res_diff_R = zeros(NS, NP, Nx)
        for nS in 1:NS, nP in 1:NP, nx in 1:Nx
        # nS = 1; nP = 3; nx = 1
            if nP >= 2 #No left for 1st paertition
                prim_res_diff_L[nS, nP, nx] = star_SP_x0[nS, nP, nx] - z_diff[nS,nP-1,nx]
            end

            if nP <= NP - 1 #No right for Last
                prim_res_diff_R[nS, nP, nx] = star_SP_x[nS, nP, nx, end] - z_diff[nS, nP, nx]
            end
        end
            
        prim_res_des
        prim_res_diff_L
        prim_res_diff_R
    #endregion
    prim_Residual_Norm = ( sum(prim_res_des[nS,nP,ndes]^2 + prim_res_diff_L[nS,nP,nx]^2 + prim_res_diff_R[nS,nP,nx]^2   for nS in 1:NS, nP in 1:NP,   ndes in 1:Ndes, nx in 1:Nx ) )^0.5

    #region ->#*Dual Residuals
        #Design
        dual_res_des = zeros(Ndes)
        for ndes in 1:Ndes
            dual_res_des[ndes]  = rho*(z_des[ndes]  - plot_z_des[end][ndes] )
        end

        #Differential
        dual_res_diff = zeros(NS, NP-1, Nx)
        for nS in 1:NS, nP in 1:NP-1, nx in 1:Nx 
            dual_res_diff[nS, nP, nx] = rho*(  z_diff[nS, nP, nx]- plot_z_diff[end][nS, nP, nx] )
        end

        dual_res_des
        dual_res_diff
    #endregion
    dual_Residual_Norm = ( sum(dual_res_des[ndes]^2       + dual_res_diff[nS,nP,nx]^2                                   for nS in 1:NS, nP in 1:NP-1, ndes in 1:Ndes, nx in 1:Nx) )^0.5

    #region - Oldcode
        # prim_res_des = NaN*zeros(1,3)
        #     prim_res_des[1,1]  = star_des1[1]  - z_des[1]
        #     prim_res_des[1,2]  = star_des2[1]  - z_des[1]
        #     prim_res_des[1,3]  = star_des3[1]  - z_des[1]
            
        
        # prim_res_diff = zeros(1,3,2)  #Dimension NS X NP X NEnd
        #     prim_res_diff[1,1,2] = star_x1[1, end, end] - z_diff[1]
        #     prim_res_diff[1,2,1] = star_x02[1]          - z_diff[1]

        #     prim_res_diff[1,2,2] = star_x2[1, end, end] - z_diff[2]
        #     prim_res_diff[1,3,1] = star_x03[1]          - z_diff[2]

        # prim_res =  (sum(prim_res_des[1,nP]^2 for nP in 1:3) + sum( prim_res_diff[1,nP,nEnd]^2 for nP in 1:3, nEnd in 1:2) )^0.5

        # dual_res_des  = rho*(z_des[1]  - plot_z_des[end])
        # dual_res_diff = rho*(z_diff[1] - plot_z_diff[end])
        # dual_res = (dual_res_des[1]^2 + dual_res_diff[1]^2)^0.5

        # μ1_des[1] = μ1_des[1] + rho*(star_des1[1] - z_des[1])
        # μ2_des[1] = μ2_des[1] + rho*(star_des2[1] - z_des[1])
        # μ3_des[1] = μ3_des[1] + rho*(star_des3[1] - z_des[1])


        # μ1_diff[1] = μ1_diff[1] + rho*(star_x1[1, end, end]   - z_diff[1])
        # μ2_diff[1] = μ2_diff[1] + rho*(star_x02[1]            - z_diff[1])
        # μ2_diff[2] = μ2_diff[2] + rho*(star_x2[1, end, end]   - z_diff[2])
        # μ3_diff[1] = μ3_diff[1] + rho*(star_x03[1]            - z_diff[2])

        
            #Storing in Plots
                
                # #scaled values
                # append!(plot_z_des, z_des)
                # append!(plot_z_diff, z_diff)
                
                # append!(plot_des1, star_des1)
                # append!(plot_x01, star_x01)
                # plot_x1 = cat(plot_x1, transpose(star_x1), dims = 2)
                # append!(plot_μ1_des, μ1_des)
                # append!(plot_μ1_diff, μ1_diff)
                # append!(plot_Obj1, star_Obj1)

                #     append!(plot_des2, star_des2)
                #     append!(plot_x02, star_x02)
                #     plot_x2 = cat(plot_x2, transpose(star_x2), dims = 2)
                #     append!(plot_μ2_des, μ2_des)
                #     append!(plot_μ2_diff, μ2_diff)
                #     append!(plot_Obj2, star_Obj2)

                # #Unscaled values
                # append!(plot_V_tes1, star_V_tes1)
                # plot_T_tes1 = cat(plot_T_tes1, star_T_tes1, dims = 2)

                #     append!(plot_V_tes2, star_V_tes2)
                #     plot_T_tes2 = cat(plot_T_tes2, star_T_tes2, dims = 2)


    #endregion

    #region ->#*μ Update
        for nS in 1:NS, nP in 1:NP, ndes in 1:Ndes
            μ_des[nS,nP,ndes] = μ_des[nS,nP,ndes] + rho*(prim_res_des[nS,nP,ndes])
        end

        for nS in 1:NS, nP in 1:NP, nx in 1:Nx
            μ_diff_L[nS,nP,nx] =  μ_diff_L[nS,nP,nx] + rho*(prim_res_diff_L[nS,nP,nx])
            μ_diff_R[nS,nP,nx] =  μ_diff_R[nS,nP,nx] + rho*(prim_res_diff_R[nS,nP,nx])
        end

    #endregion

    #region ->#* rho update heuristic
        if prim_Residual_Norm       > 10*dual_Residual_Norm
            rho = rho*2
        elseif dual_Residual_Norm   > 10*prim_Residual_Norm
            rho = rho/2
        else
        end
    #endregion

    #region ->#* Appending for Plotting

        push!(plot_rho,       copy(rho))
        push!(plot_z_des,     copy(z_des) )
        push!(plot_z_diff,    copy(z_diff))

        push!(plot_Obj,       star_SP_Obj)
        push!(plot_des,       star_SP_des)
        push!(plot_x,         star_SP_x)
        push!(plot_x0,        star_SP_x0)

        push!(plot_μ_des,     copy(μ_des))
        push!(plot_μ_diff_L,  copy(μ_diff_L))
        push!(plot_μ_diff_R,  copy(μ_diff_R))

        push!(plot_V_tes,     star_SP_V_tes)
        push!(plot_T_tes,     star_SP_T_tes)
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

plot_x0
plot_z_des
plot_z_diff

plot_μ_des
plot_μ_diff_L   #!check if adjacent ones sum to zero (with SP3)
plot_μ_diff_R

##* Calculating ADMM - Summary Stats

    #region ->Primal Residual
        plot_prim_res_des = []
        plot_prim_res_diff_L = []
        plot_prim_res_diff_R = []
        plot_prim_Residual_Norm = []
        for nIter in 1:NIter
            # nIter = 1

            #Design
            prim_res_des    = zeros(NS, NP, Ndes)
            for nS in 1:NS, nP in 1:NP, ndes in 1:Ndes
                prim_res_des[nS, nP, ndes] = plot_des[nIter][nS, nP, ndes] - plot_z_des[nIter][ndes]
            end
            
            #Differential
            prim_res_diff_L = zeros(NS, NP, Nx)
            prim_res_diff_R = zeros(NS, NP, Nx)
            for nS in 1:NS, nP in 1:NP, nx in 1:Nx
                # nS = 1; nP = 3; nx = 1
                if nP >= 2 #No left for 1st paertition
                    prim_res_diff_L[nS, nP, nx] = plot_x0[nIter][nS, nP, nx]     - plot_z_diff[nIter][nS, nP-1, nx]
                    
                end

                if nP <= NP - 1 #No right for Last
                    prim_res_diff_R[nS, nP, nx] = plot_x[nIter][nS, nP, nx, end] - plot_z_diff[nIter][nS, nP, nx]
                end
            end

            #2-Norm of primal residual
            prim_Residual_Norm = (sum( prim_res_des[nS,nP,ndes]^2 + prim_res_diff_L[nS,nP,nx]^2 + prim_res_diff_R[nS,nP,nx]^2    for nS in 1:NS, nP in 1:NP, ndes in 1:Ndes, nx in 1:Nx) )^0.5
            
            push!(plot_prim_res_des,        prim_res_des)
            push!(plot_prim_res_diff_L,     prim_res_diff_L)
            push!(plot_prim_res_diff_R,     prim_res_diff_R)  
            push!(plot_prim_Residual_Norm,  prim_Residual_Norm)
                        
        end

        plot_prim_res_des
        plot_prim_res_diff_L
        plot_prim_res_diff_R
        plot_prim_Residual_Norm

    #endregion

    #region ->Dual Residual
        plot_dual_res_des  = []
        plot_dual_res_diff = []
        plot_dual_Residual_Norm = []
        for nIter in 1:NIter
            # nIter = 2

            dual_res_des = zeros(Ndes)
            for ndes in 1:Ndes
                dual_res_des[ndes] = plot_rho[nIter]* (plot_z_des[nIter+1][ndes] - plot_z_des[nIter][ndes] )
            end

            dual_res_diff = zeros(NS,NP-1,Nx)*NaN
            for nS in 1:NS, nP in 1:NP-1, nx in 1:Nx
                dual_res_diff[nS,nP,nx] =  plot_rho[nIter]* (plot_z_diff[nIter+1][nS,nP,nx] - plot_z_diff[nIter][nS,nP,nx])
            end
            
            #2-Norm of residual
            dual_Residual_Norm = (sum(dual_res_des[ndes]^2 + dual_res_diff[nS,nP,nx]^2    for nS in 1:NS, nP in 1:NP-1, ndes in 1:Ndes, nx in 1:Nx)  )^0.5

            push!(plot_dual_res_des,    dual_res_des)
            push!(plot_dual_res_diff,   dual_res_diff)
            push!(plot_dual_Residual_Norm,   dual_Residual_Norm)

        end
                
        plot_dual_res_des
        plot_dual_res_diff
        plot_dual_Residual_Norm

    #endregion

    #region ->Augmented Terms
        plot_Aug_term   = []
        for nIter in 1:NIter   
            # nIter = 2

            Aug_term      = zeros(NS,NP)*NaN
            Aug_term_des  = zeros(NS,NP)*NaN
            Aug_term_diff = zeros(NS,NP)*NaN
            for nS in 1:NS, nP in 1:NP            
                # nS = 1; nP = 2

                #Design Augmented terms
                Aug_term_des[nS,nP]      = sum(plot_μ_des[nIter][nS,nP,ndes]*(plot_prim_res_des[nIter][nS,nP,ndes])            + rho*(plot_prim_res_des[nIter][nS,nP,ndes])^2            for ndes in 1:Ndes)

                #differential augmented terms -> according to leftmost, Rightmost or Inbetween
                if nP == 1      #*only right

                    Aug_term_diff[nS,nP] = sum(plot_μ_diff_R[nIter][nS,nP,nx]*(plot_prim_res_diff_R[nIter][nS, nP, nx])       + rho*(plot_prim_res_diff_R[nIter][nS, nP, nx])^2          for nx in 1:Nx)

                elseif nP == NP #*Only Left
                    
                    Aug_term_diff[nS,nP] = sum(plot_μ_diff_L[nIter][nS,nP,nx]*(plot_prim_res_diff_L[nIter][nS,nP,nx])         + rho*(plot_prim_res_diff_L[nIter][nS,nP,nx])^2            for nx in 1:Nx)

                else            #*Right and Left

                    Aug_term_diff[nS,nP] = sum(plot_μ_diff_R[nIter][nS,nP,nx]*(plot_prim_res_diff_R[nIter][nS, nP, nx])       + rho*(plot_prim_res_diff_R[nIter][nS, nP, nx])^2          for nx in 1:Nx) +
                                           sum(plot_μ_diff_L[nIter][nS,nP,nx]*(plot_prim_res_diff_L[nIter][nS,nP,nx])         + rho*(plot_prim_res_diff_L[nIter][nS,nP,nx])^2            for nx in 1:Nx)

                end
            
                Aug_term[nS,nP] = Aug_term_des[nS,nP] + Aug_term_diff[nS,nP]
            
            end

            push!(plot_Aug_term, Aug_term)                 

        end 

        plot_Aug_term

                #These must be undefined (because of local scopes)
                # Aug_term
                # Aug_term_des
                # Aug_term_diff

    #endregion

    #region ->Original Objective function value
        plot_star_f    = []
        for nIter in 1:NIter  
            # nIter = 1 
            
            star_f = zeros(NS,NP)*NaN
            for nS in 1:NS, nP in 1:NP
                # nS = 1; nP = 1

                star_f[nS,nP] = plot_Obj[nIter][nS,nP]     -   plot_Aug_term[nIter][nS,nP] 
            end
            
            push!(plot_star_f, star_f)
        end

        plot_star_f

                #These must be undefined (because of local scopes)
                # star_f

    #endregion

    plot_prim_res_des
    plot_prim_res_diff_L
    plot_prim_res_diff_R
    plot_prim_Residual_Norm

    plot_dual_res_des
    plot_dual_res_diff
    plot_dual_Residual_Norm

    plot_Aug_term
    
    plot_star_f

##* Plot ADMM itreations Summary

    #choose backend for plots
    plotly()
    # gr()

    #Profiles - Differential States
    nS1, nP1 = 1,1
    nS2, nP2 = 1,2

        p11 = plot( collect(0:30),  vcat(plot_x0[end][nS1,nP1,1]*100, plot_T_tes[end][nS1,nP1,:]) ,          label = "Ttes - SP1", ylim = [55,70])
        p11 = plot!(collect(30:60), vcat(plot_x0[end][nS2,nP2,1]*100, plot_T_tes[end][nS2,nP2,:]) ,          label = "Ttes - SP2", title = "Diff state, rho = $rho", ylim = [55,70])

    #Iteration Profiles - Design Variable
    plot_Iter = collect(1:NIter)
        p21 = plot( plot_Iter, [plot_V_tes[nIter][1]  for nIter in 1:NIter],                    label = "des - SP1")
        p21 = plot!(plot_Iter, [plot_V_tes[nIter][2]  for nIter in 1:NIter],                    label = "des - SP2")
        p21 = plot!(plot_Iter, star_V_tes.*ones(NIter),                                         label = "des Centr",    title = "Des var, rho = $rho" )

    #Primal and Dual Infeasibility
        p31 = plot( plot_Iter, [plot_prim_Residual_Norm[nIter]  for nIter in 1:NIter],          label  = "primal infeasibility", yscale = :log10)
        p31 = plot!(plot_Iter, [plot_dual_Residual_Norm[nIter]  for nIter in 1:NIter],          label = "dual infeasibility",    yscale = :log10, title = "Infeasibility, rho = $rho")

    #Objective Function Values
        p41 = plot(plot_Iter,  [ sum(plot_star_f[nIter][nS,nP]  for nS in 1:NS, nP in 1:NP) for nIter in 1:NIter] , label = "f - Distribuited")
        p41 = plot!(plot_Iter, star_Obj.*ones(NIter),                                                               label = "f - Central Sol'n" )
        p41 = plot!(plot_Iter, [ sum(plot_Obj[nIter][nS,nP]     for nS in 1:NS, nP in 1:NP) for nIter in 1:NIter],  label = "Augmented Lagrange", title = "Objective, rho = $rho")
    
        # p51 = plot(plot_Iter, plot_f1_star,                 label ="f1 star")
        # p51 = plot!(plot_Iter, plot_f2_star,                label ="f2 star")
        # p51 = plot!(plot_Iter, plot_Aug_terms1,             label = "f1 Aug Terms")
        # p51 = plot!(plot_Iter, plot_Aug_terms2,             label = "f2 Aug Terms", title = "SP1 objective, rho = $rho")

        p61 = plot( plot_Iter, [plot_μ_des[nIter][1,1]      for nIter in 1:NIter],              label = "μ1_des")
        p61 = plot!(plot_Iter, [plot_μ_des[nIter][1,2]      for nIter in 1:NIter],              label = "μ2_des")
        p61 = plot!(plot_Iter, [plot_μ_diff_R[nIter][1,1]   for nIter in 1:NIter],              label = "μ1 diff_R")
        p61 = plot!(plot_Iter, [plot_μ_diff_L[nIter][1,2]   for nIter in 1:NIter],              label = "μ2 diff_L",    title = "multipliers, rho = $rho")


        p71 = plot(plot_Iter, plot_rho[2:end],              title = "Rho updates")


##* Display Plots
p11
p21
p31
p41
# p51
p61
display(Tot_time_in_ADMM)


