
using Plots
using JuMP
using Ipopt

##todo - changing q_dh in file does not change in this file. Reloading module how?
include("parameters.jl")
      using Main.Bounds
                        # :      x0,      z0,     u0,
                        #       Nx,      Nz,     Nu, 
                        #       ls_x,    ls_z,   ls_u, 
                        #       us_x,    us_z,   us_u  , x0_us  #scaled variables 

      using Main.Model_parameters

function Collocation_Matrix()
      #Radau
      t1 = 0.155051
      t2 = 0.644949
      t3 = 1.0

      M1 = [
            t1 1 / 2 * t1^2 1 / 3 * t1^3
            t2 1 / 2 * t2^2 1 / 3 * t2^3
            t3 1 / 2 * t3^2 1 / 3 * t3^3
            ]
      M2 = [
            1 t1 t1^2
            1 t2 t2^2
            1 t3 t3^2
            ]

      M = M1 * inv(M2)
      return M
end

# function Solve_OCP(x0, Tf, Is_ISS)

            #region -> Value of Arguments for Debugging
            # x0 = [60.0]
            Tf = 30.0
            Is_SS = 0          #If want to find Steady State set = 1

            #endregion

                  #region -> Setting Initial guesses and Dimensions
                  # z0 = [60.0; 60.0; 60.0]
                  # u0 = [0.0; 0.0]

                  dx0_us  = 0 * x0
                  alg0_us = 0 * z0
                  # Nx = size(x0, 1)
                  # Nz = size(z0, 1)
                  # Nu = size(u0, 1)
                  #endregion



      ##*OCP Parameters
      T0 = 0.0
      NFE = 30
      dt = (Tf - T0) / NFE
      NCP = 3


##* Defining Solver
      model1 = Model(with_optimizer(Ipopt.Optimizer))

      #region -> Variables and Objective
            ## Declare Variables
            @variable(model1, x[1:Nx, 1:NFE, 1:NCP])
            @variable(model1, dx_us[1:Nx, 1:NFE, 1:NCP])

            @variable(model1, z[1:Nz, 1:NFE, 1:NCP])
            # @variable(model1, alg_us[1:Nz, 1:NFE, 1:NCP])

            @variable(model1, u[1:Nu, 1:NFE])

                  #region-> #? Set Variable Bounds AND Initial Guesses (scaled)
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP     
                        set_lower_bound(x[nx, nfe, ncp], 0)
                        # set_upper_bound(x[nx, nfe, ncp], 999)

                        set_lower_bound(z[nz, nfe, ncp], 0)
                        #set_upper_bound(z[nz, nfe, ncp], 999)

                        #set_lower_bound(dx_us[nx, nfe, ncp], 0)
                        #set_upper_bound(dx_us[nx, nfe, ncp], 999)

                        #set_lower_bound(alg_us[nz, nfe, ncp], 0)
                        #set_upper_bound(alg_us[nz, nfe, ncp], 999)

                        set_lower_bound(u[nu, nfe], 0)
                        set_upper_bound(u[1,  nfe], 1)
                  end

                        #Initial Guesses (Scaled)
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP
                        set_start_value(x[nx, nfe, ncp],          x0[nx])
                        set_start_value(z[nz, nfe, ncp],          z0[nz])
                        set_start_value(dx_us[nx, nfe, ncp],      dx0_us[nx])
                        # set_start_value(alg_us[nz, nfe, ncp],     alg0_us[nz])
                        set_start_value(u[nu, nfe],               u0[nu])


                  end
                  #endregion

                  #region-> Expressions for Unscaling Variables (easier to write DAE Equation) #?add variables and change indices 
                  @NLexpressions(model1, begin

                        T_tes[nfe in 1:NFE, ncp in 1:NCP],  x[1, nfe, ncp]            *  (us_x[1] - ls_x[1]) + ls_x[1]

                        T_b[nfe in 1:NFE, ncp in 1:NCP],    z[1, nfe, ncp]            *  (us_z[1] - ls_z[1]) + ls_z[1]
                        T_phb[nfe in 1:NFE, ncp in 1:NCP],  z[2, nfe, ncp]            *  (us_z[2] - ls_z[2]) + ls_z[2]
                        T_whb[nfe in 1:NFE, ncp in 1:NCP],  z[3, nfe, ncp]            *  (us_z[3] - ls_z[3]) + ls_z[3]

                        α[nfe in 1:NFE],        u[1, nfe]                             *  (us_u[1] - ls_u[1]) + ls_u[1] 
                        Q_phb[nfe in 1:NFE],    u[2, nfe]                             *  (us_u[2] - ls_u[2]) + ls_u[2] 

                  end)
                  #endregion

            ## Objective
            @NLobjective(model1, Min, sum( Q_phb[nfe] for nfe in 1:NFE ) )    #Duty in KJ

      #endregion

##* Defining Constraints
     #region -> Constraints

           @NLconstraints(model1, begin
            #?Defining the model ODEs in each line
            Constr_ODE1[nfe in 1:NFE, ncp in 1:NCP], dx_us[1, nfe, ncp]      == α[nfe]*q_dh*3600*( T_whb[nfe,ncp] - T_tes[nfe,ncp] )/V_tes            #DegC/hr
            #In case of more states - pattern
            #Constr_ODE999[nfe=1:NFE, ncp=1:NCP], dx[999,nfe,ncp] ==
            end)

            @NLconstraints(model1, begin
            #?Defining Model Algebraic Equations in each line
            Constr_Alg1[nfe in 1:NFE, ncp in 1:NCP], T_b[nfe, ncp]      == α[nfe]*T_tes[nfe,ncp]      + (1-α[nfe])*T_whb[nfe,ncp]
            Constr_Alg2[nfe in 1:NFE, ncp in 1:NCP], T_phb[nfe, ncp]    == T_b[nfe,ncp]               + Q_phb[nfe]/( q_dh*ρ_dh*Cp_dh )
            Constr_Alg3[nfe in 1:NFE, ncp in 1:NCP], T_whb[nfe, ncp]    == T_dh_ret                   + Q_whb[nfe]/( q_dh*ρ_dh*Cp_dh) 
            #In case of more states - pattern
            #Constr_Alg999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
            end)


            @NLconstraints(model1, begin
                  #?Defining any Inequality Constraints in each line
                  Constr_Ineq1[nfe in 1:NFE, ncp in 1:NCP], T_phb[nfe, ncp] >= T_dh_minSup
                  #In case of more states - pattern
                  #Constr_Ineq999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
            end)

                  #region-> std code -> Collocation Equations
                        collMat = Collocation_Matrix()
                        @NLconstraints(model1, begin
                              #Collocation Equation for Differential Equations (scaled form)
                              #t = 0
                              Constr_Coll_Diff0[nx in 1:Nx,  nfe = 1, ncp in 1:NCP],        x[nx, nfe, ncp]    == x0[nx]               + dt * sum(collMat[ncp, i] * dx_us[nx, nfe, i] for i = 1:NCP) /(us_x[nx] - ls_x[nx])
                              #t = 1 ... (N-1)
                              Constr_Coll_Diff[nx in 1:Nx,   nfe in 2:NFE, ncp = 1:NCP],    x[nx, nfe, ncp]    == x[nx, nfe-1, NCP]    + dt * sum(collMat[ncp, i] * dx_us[nx, nfe, i] for i = 1:NCP) /(us_x[nx] - ls_x[nx])
                        end)
                        
                  #endregion-> std code     

                  
      #endregion

##* Solve the model
      optimize!(model1)
      JuMP.termination_status(model1)
      JuMP.solve_time(model1::Model)

      #region-> Generic accessing
      # star_x = JuMP.value.(x[:, :, NCP])
      # star_x = cat(x0, star_x, dims = 2)
      # star_u = JuMP.value.(u[:, :])
      # star_z = JuMP.value.(z[:, :, NCP])
      star_dx_us = JuMP.value.(dx_us[:,:,NCP])
      #endregion
     
      #Getting Varables directly from expressions
      star_T_tes = JuMP.value.(T_tes[:, NCP])
                  star_T_tes = cat(x0[1]*100, star_T_tes, dims = 1)
      
      star_T_b   = JuMP.value.(T_b[:, NCP])
      star_T_phb   = JuMP.value.(T_phb[:, NCP])
      star_T_whb   = JuMP.value.(T_whb[:, NCP])
      
      star_α      = JuMP.value.(α[:])
      star_Q_phb = JuMP.value.(Q_phb[:])
     

##* Plot Solution
      #region -> Plotting Solution

            t_plot = collect(T0:dt:Tf)    #Returns NFE+1 dimensional vector
            
            #choose backend for plots
            # plotly()
            gr()

                  #Differential States
                  p11 = plot(t_plot, star_T_tes,                  label = "T_tes",
                                                                                    marker = true)
                  #Manipulated Variables
                  p12 = plot(t_plot[1:end-1], star_α,             label = "α",
                                                                                    marker = true, linetype = :steppost)
                  p13 = plot(t_plot[1:end-1], star_Q_phb,        label = "Q_phb",
                                                                                    marker = true, linetype = :steppost)
                  #Algebraic States
                  p14 = plot(t_plot[2:end], star_T_b,             label = "T_b", 
                                                                                    marker = true)
                  p14 = plot!(t_plot[2:end], star_T_phb,          label = "T_phb", 
                                                                                    marker = true)
                  p14 = plot!(t_plot[2:end], star_T_whb,          label = "T_whb", 
                                                                                    marker = true)

            fig1 = plot(p11, p12, p13, p14, layout = (4, 1))

      #                         #Differential States
      #                         p11 = plot(t_plot, star_x[1, :],                label = "T_tes",
      #                         marker = true)
      #                         #Manipulated Variables
      #                         p12 = plot(t_plot[1:end-1], star_u[1, :],       label = "α",
      #                                                       marker = true, linetype = :steppost)
      #                         p12 = plot!(t_plot[1:end-1], star_u[2, :],      label = "Q_phb",
      #                                                       marker = true, linetype = :steppost)
      #                         #Algebraic States
      #                         p13 = plot(t_plot[2:end], star_z[1,:],          label = "T_b", 
      #                                                       marker = true)
      #                         p13 = plot!(t_plot[2:end], star_z[2,:],         label = "T_phb", 
      #                                                       marker = true)
      #                         p13 = plot!(t_plot[2:end], star_z[3,:],         label = "T_whb", 
      #                                                       marker = true)

      #                         fig1 = plot(p11, p12, p13, layout = (3, 1))
      # #endregion


#       return star_u[:, 1]


# end


