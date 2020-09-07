#=
28/8/20 - Exercise 5 from MPC modeule 2019 version at NTNU.
Creates collocation matrix of 3rd order Radau polynomial
Solves OCP taking in initiial differential state, previous MV (for penalizing input moves) and setpoint trajectory
=#

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

function Solve_OCP(x0, u₋₁, y_sp)

            #region -> Value of Arguments for Debugging
            # x0 = [1.0; 1.0]
            # u₋₁ = 0.3
            # y_sp = append!(1.5302*ones(10), 0.9951*ones(10) )

            #endregion

      ##*Model Parameters 
      x2_f = 4
      km = 0.12
      k1 = 0.4545
      Y = 0.4
      μ_max = 0.53
      D_min = 0
      D_max = 1

      ##*OCP Parameters
      T0 = 0.0
      Tf = 30.0
      NFE = 30
      dt = (Tf - T0) / NFE
      NCP = 3

                  #region -> Setting Initial guesses and Dimensions
                  q0 = 0
                  z0 = [0.0]
                  u0 = [0.0]

                  dx0 = 0 * x0
                  alg0 = 0 * z0
                  dq0 = 0
                  Nx = size(x0, 1)
                  Nz = size(z0, 1)
                  Nu = size(u0, 1)
                  #endregion

      ##* Defining Solver
      model1 = Model(with_optimizer(Ipopt.Optimizer))

            #region -> Variables and Objective
            ## Declare Variables
            @variable(model1, x[1:Nx, 1:NFE, 1:NCP])
            @variable(model1, dx[1:Nx, 1:NFE, 1:NCP])

            @variable(model1, q[1, 1:NFE, 1:NCP])
            @variable(model1, dq[1, 1:NFE, 1:NCP])

            @variable(model1, z[1:Nz, 1:NFE, 1:NCP])
            @variable(model1, alg[1:Nz, 1:NFE, 1:NCP])

            @variable(model1, u[1:Nu, 1:NFE])

                  #region -> Set Variable Bounds AND Initial Guesses 
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP     
                        set_lower_bound(x[nx, nfe, ncp], 0)
                        set_upper_bound(x[1, nfe, ncp], 4.5)

                        set_lower_bound(z[nz, nfe, ncp], 0)
                        #set_upper_bound(z[nz, nfe, ncp], 999)

                        #set_lower_bound(dx[nx, nfe, ncp], 0)
                        #set_upper_bound(dx[nx, nfe, ncp], 999)

                        #set_lower_bound(alg_coll[nz, nfe, ncp], 0)
                        #set_upper_bound(alg_coll[nz, nfe, ncp], 999)

                        set_lower_bound(u[nu, nfe], D_min)
                        set_upper_bound(u[nu, nfe], D_max)
                  end

                  #Set Initial Guesses
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP
                        set_start_value(x[nx, nfe, ncp], x0[nx])
                        set_start_value(z[nz, nfe, ncp], z0[nz])
                        set_start_value(dx[nx, nfe, ncp], dx0[nx])
                        set_start_value(alg[nz, nfe, ncp], alg0[nz])
                        set_start_value(u[nu, nfe], u0[nu])
                        set_start_value(q[1, nfe, ncp], q0)
                        set_start_value(dq[1, nfe, ncp], dq0)

                  #endregion
                  end

            ## Objective
            @NLobjective(model1, Min, q[1, NFE, NCP])
            
            #endregion


      #region -> Constraints

      @NLconstraints(model1, begin
            #Defining the model ODEs in each line
            Constr_ODE1[nfe in 1:NFE, ncp in 1:NCP], dx[1, nfe, ncp]      == x[1, nfe, ncp] * (z[1, nfe, ncp] - u[1, nfe])
            Constr_ODE2[nfe in 1:NFE, ncp in 1:NCP], dx[2, nfe, ncp]      == u[1, nfe]      * (x2_f - x[2, nfe, ncp])           - z[1, nfe, ncp] * x[1, nfe, ncp] / Y
            #In case of more states - pattern
            #Constr_ODE999[nfe=1:NFE, ncp=1:NCP], dx[999,nfe,ncp] ==
      end)

      @NLconstraints(model1, begin
            #Collocation for Quadrature/ Objective Function
            #t = 0
            Constr_quad_dot0[nfe = 1    , ncp in 1:NCP], dq[1, nfe, ncp]  == (x[1,nfe,ncp] - y_sp[nfe])^2 + 0.5*(u[1,nfe]  - u₋₁[1]    )^2
            #t = 1..N
            Constr_quad_dot[nfe in 2:NFE, ncp in 1:NCP], dq[1, nfe, ncp]  == (x[1,nfe,ncp] - y_sp[nfe])^2 + 0.5*(u[1,nfe] - u[1,nfe-1])^2
            #In case of objective depends on end state only -> dq = dx[2] etc
      end)


      @NLconstraints(model1, begin
            #Defining Model Algebraic Equations in each line
            Constr_Alg1[nfe in 1:NFE, ncp in 1:NCP], z[1, nfe, ncp] == (  μ_max * x[2, nfe, ncp] / (km + x[2, nfe, ncp] + k1 * x[2, nfe, ncp]^2) )
            #In case of more states - pattern
            #Constr_Alg999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
      end)

      @NLconstraints(model1, begin
            #Defining Inequality Constraints in each line
            Constr_Ineq1[nfe in 1:1   ], u[1,nfe] - u₋₁[1]       <= 0.05
            Constr_Ineq2[nfe in 2:NFE ], u[1,nfe] - u[1,nfe-1]   <= 0.05
            #In case of more states - pattern
            #Constr_Ineq999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
      end)


            #region-> std code -> Collocation Equation
            ## Creating a Radau collocation Matrix for NCP = 3
            collMat = Collocation_Matrix()
            @NLconstraints(model1, begin
                  #Collocation Equation for Differential Equations AND Objective Function
                  #t = 0
                  Constr_Coll_Diff0[nx in 1:Nx,  nfe = 1, ncp in 1:NCP],        x[nx, nfe, ncp]    == x0[nx]              + dt * sum(collMat[ncp, i] * dx[nx, nfe, i] for i = 1:NCP)
                  Constr_Coll_quad0[             nfe = 1, ncp in 1:NCP],        q[1, nfe, ncp]     == q0                  + dt * sum(collMat[ncp, i] * dq[1, nfe, i]  for i = 1:NCP)
                  #t = 1 ... (N-1)
                  Constr_Coll_Diff[nx in 1:Nx,   nfe in 2:NFE, ncp = 1:NCP],    x[nx, nfe, ncp]     == x[nx, nfe-1, NCP]   + dt * sum(collMat[ncp, i] * dx[nx, nfe, i] for i = 1:NCP)
                  Constr_Coll_quad[              nfe in 2:NFE, ncp = 1:NCP],    q[1, nfe, ncp]      == q[1, nfe-1, NCP]    + dt * sum(collMat[ncp, i] * dq[1, nfe, i]  for i = 1:NCP)
            end)
            #endregion std code

      #endregion

      ## Solve the model
      optimize!(model1)
      JuMP.termination_status(model1)
      JuMP.solve_time(model1::Model)

      star_x = JuMP.value.(x[:, :, NCP])
      star_x = cat(x0, star_x, dims = 2)
      star_u = JuMP.value.(u[:, :])
      star_z = JuMP.value.(z[:, :, NCP])

      #region -> Plotting Solution

            # t_plot = collect(T0:dt:Tf)    #Returns NFE+1 dimensional vector
            #
            # #choose backend for plots
            # plotly()
            #
            # p11 = plot(t_plot, star_x[1, :], label = "x1")
            # p11 = plot!(t_plot, star_x[2, :], label = "x2")
            # p11 = plot!(t_plot[1:end-1], y_sp[:], label = "y_sp",
            #                                                             linetype = :steppost, linestyle = :dash)
            # p12 = plot(t_plot[1:end-1], star_u[1, :], label = "u1",
            #                                                             linetype = :steppost, marker = false)
            # fig1 = plot(p11, p12, layout = (2, 1))

      #endregion


      return star_u[:, 1]


end


