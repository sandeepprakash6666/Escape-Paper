
using Plots
using JuMP
using Ipopt

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
            x0 = [60.0]
            Tf = 30.0
            Is_SS = 0          #If want to find Steady State set = 1

            #endregion

                  #region -> Setting Initial guesses and Dimensions
                  z0 = [60.0; 60.0; 60.0]
                  u0 = [0.0; 0.0]

                  dx0 = 0  * x0
                  alg0 = 0 * z0
                  Nx = size(x0, 1)
                  Nz = size(z0, 1)
                  Nu = size(u0, 1)
                  #endregion

      ##*Model Parameters 
      ρ_dh = 1000
      Cp_dh = 4.18
      q_dh = 1

      ρ_wh = 1000
      Cp_wh = 4.18
      q_wh = 1

      T_dh_ret = 60
      T_dh_minSup = 60

      Q_whb = 0
      V_tes = 500

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
            @variable(model1, dx[1:Nx, 1:NFE, 1:NCP])

            @variable(model1, z[1:Nz, 1:NFE, 1:NCP])
            @variable(model1, alg[1:Nz, 1:NFE, 1:NCP])

            @variable(model1, u[1:Nu, 1:NFE])

                  #region -> Set Variable Bounds AND Initial Guesses (scaled)
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP     
                        set_lower_bound(x[nx, nfe, ncp], 0)
                        # set_upper_bound(x[nx, nfe, ncp], 999)

                        set_lower_bound(z[nz, nfe, ncp], 0)
                        #set_upper_bound(z[nz, nfe, ncp], 999)

                        #set_lower_bound(dx[nx, nfe, ncp], 0)
                        #set_upper_bound(dx[nx, nfe, ncp], 999)

                        #set_lower_bound(alg[nz, nfe, ncp], 0)
                        #set_upper_bound(alg[nz, nfe, ncp], 999)

                        set_lower_bound(u[1, nfe], 0)
                        set_upper_bound(u[1, nfe], 1)
                        set_lower_bound(u[2, nfe], 0)
                        # set_upper_bound(u[2, nfe], 0)                       
                  end

                  #todo - Set Initial Guesses as scaled
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP
                        set_start_value(x[nx, nfe, ncp],    x0[nx])
                        set_start_value(z[nz, nfe, ncp],    z0[nz])
                        set_start_value(dx[nx, nfe, ncp],   dx0[nx])
                        set_start_value(alg[nz, nfe, ncp],  alg0[nz])
                        set_start_value(u[nu, nfe],         u0[nu])


                  #endregion
                  end

                  #todo Make expressions for scaling and naming variables (for ease of writing eqns)
                  #region-> Expressions for Unscaling Variables
                  # @NLexpressions(model1, begin
                  
                        # Z1[nfe in 1:NFE, ncp in 1:NCP], 2.0*z[1, nfe, ncp]

                  # end)
                  #endregion


            ## Objective
            @NLobjective(model1, Min, sum( u[2,nfe] for nfe in 1:NFE ) )

            #endregion

##* Defining Constraints
     #region -> Constraints

      @NLconstraints(model1, begin
            #Defining the model ODEs in each line
            Constr_ODE1[nfe in 1:NFE, ncp in 1:NCP], dx[1, nfe, ncp]      == u[1, nfe]*q_dh*(z[3,nfe,ncp] - x[1,nfe,ncp] )/V_tes
            #In case of more states - pattern
            #Constr_ODE999[nfe=1:NFE, ncp=1:NCP], dx[999,nfe,ncp] ==
      end)

      @NLconstraints(model1, begin
            #Defining Model Algebraic Equations in each line
            Constr_Alg1[nfe in 1:NFE, ncp in 1:NCP], z[1, nfe, ncp] == u[1,nfe]*x[1,nfe,ncp]    + (1-u[1,nfe])*z[3,nfe,ncp]
            Constr_Alg2[nfe in 1:NFE, ncp in 1:NCP], z[2, nfe, ncp] == z[1,nfe,ncp]             + u[2,nfe]/( q_dh*ρ_dh*Cp_dh )
            Constr_Alg3[nfe in 1:NFE, ncp in 1:NCP], z[3, nfe, ncp] == T_dh_ret                 + Q_whb/   ( q_dh*ρ_dh*Cp_dh) 
            #In case of more states - pattern
            #Constr_Alg999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
      end)

      @NLconstraints(model1, begin
            #Defining any Inequality Constraints in each line
            Constr_Ineq1[nfe in 1:NFE, ncp in 1:NCP], z[2, nfe, ncp] >= T_dh_minSup
            #In case of more states - pattern
            #Constr_Ineq999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
      end)

            #region-> std code -> Collocation Equations
                  ## Creating a Radau collocation Matrix for NCP = 3
                  collMat = Collocation_Matrix()
                  @NLconstraints(model1, begin
                        #Collocation Equation for Differential Equations
                        #t = 0
                        Constr_Coll_Diff0[nx in 1:Nx,  nfe = 1, ncp in 1:NCP],        x[nx, nfe, ncp]    == x0[nx]               + dt * sum(collMat[ncp, i] * dx[nx, nfe, i] for i = 1:NCP)
                        #t = 1 ... (N-1)
                        Constr_Coll_Diff[nx in 1:Nx,   nfe in 2:NFE, ncp = 1:NCP],    x[nx, nfe, ncp]     == x[nx, nfe-1, NCP]   + dt * sum(collMat[ncp, i] * dx[nx, nfe, i] for i = 1:NCP)
                  end)
                  
            #endregion std code     

            
      #endregion

##* Solve the model
      optimize!(model1)
      JuMP.termination_status(model1)
      JuMP.solve_time(model1::Model)

      star_x = JuMP.value.(x[:, :, NCP])
      star_x = cat(x0, star_x, dims = 2)
      star_u = JuMP.value.(u[:, :])
      star_z = JuMP.value.(z[:, :, NCP])

      # A = JuMP.value.(z[:, :, :])
      # B = A[1,:,:]
      # B[1,1] = 100 
      # T1[:,1] .= 100
      # T1
      # star_test

      star_dx = JuMP.value.(dx[:,:,NCP])

##* Plot Solution
      #region -> Plotting Solution

            t_plot = collect(T0:dt:Tf)    #Returns NFE+1 dimensional vector
            
            #choose backend for plots
            # plotly()
            gr()

                  #Differential States
                  p11 = plot(t_plot, star_x[1, :],                label = "T_tes",
                                                                                    marker = true)
                  #Manipulated Variables
                  p12 = plot(t_plot[1:end-1], star_u[1, :],       label = "α",
                                                                                    marker = true, linetype = :steppost)
                  p12 = plot!(t_plot[1:end-1], star_u[2, :],      label = "Q_phb",
                                                                                    marker = true, linetype = :steppost)
                  #Algebraic States
                  p13 = plot(t_plot[2:end], star_z[1,:],          label = "T_b", 
                                                                                    marker = true)
                  p13 = plot!(t_plot[2:end], star_z[2,:],         label = "T_phb", 
                                                                                    marker = true)
                  p13 = plot!(t_plot[2:end], star_z[3,:],         label = "T_whb", 
                                                                                    marker = true)

            fig1 = plot(p11, p12, p13, layout = (3, 1))

      #endregion


#       return star_u[:, 1]


# end



##*Testing some mock codes

#when indexing parts of array, link is broken
# A = [1 2 3 4 5 6]
# B = A[1:3]

# B[1]

# B[1] = 10

# A


# C = ones(3,3,3)
# D = C[1,:,:]
# D[1,1] = 100