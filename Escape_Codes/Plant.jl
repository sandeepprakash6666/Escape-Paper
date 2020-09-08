
using DifferentialEquations
using Sundials
# using NLsolve
using Plots


#region -> Initial Points and Guesses
x0 = [1.0; 0.0]
q0 = [0.0]
z0 = [0.0]
u0 = [0.1; 0.0]
tspan = (0.0, 10.0)
#endregion


##* Model Equations for Plant
function Plant_Eqns(x,z,u)
    dx = 0.0*x    #this creates a copy, not a pointer luckily
    dq = zeros(1)
    alg = 0.0*z

    #region-> Assigning names to variables for ease of writing Equations
    # T_tes = x[1]

    # T_b = z[1]
    # T_phb = z[2]
    # T_whb = z[3]

    # alpha = u[1]
    # Q_phb = u[2]
    #endregion

      #region - Plant Parameters
      # ρ_dh = 1000
      # Cp_dh = 4.18
      # q_dh = 1

      # ρ_wh = 1000
      # Cp_wh = 4.18
      # q_wh = 1

      # T_dh_ret = 30

      # Q_whb = 1000
      # V_tes = 500
      #endregion

    # ##ODEs -
    # dx = ( alpha*q_dh*(T_whb - T_tes)/V_tes      )

    # ##Objective function - term inside integration
    # # dq = u[1]*x[1]

    # ##Alg Equations
    # alg[1] = T_b    - ( alpha*T_tes + (1-alpha)*T_whb )
    # alg[2] = T_phb  - ( T_b + Q_phb/( q_dh*ρ_dh*Cp_dh)   )
    # alg[3] = T_whb  - ( T_dh_ret + Q_whb/ (q_dh*ρ_dh*Cp_dh))

      #region -> testing with Lorentz function as DAE
        ##ODEs -
        dx[1] = -0.04*x[1] + 10000*x[2]*z[1]
        dx[2] = 0.04*x[1] - 10000*x[2]*z[1] - 30000000*x[2]^2

        ##Objective function - term inside integration
        dq[1] = 0
    
        ##Alg Equations
        alg[1] = 1 - (x[1]*x[2] + z[1])
      #endregion

    return dx, dq, alg
end

#Test Plant_Eqns
# Plant_Eqns(x0, z0, u0)


##* Function to solve a DAE using SunDials integrator
function Integrate_Plant_DAE(x0, q0, z0, u0, tspan)

    #region-> Setting Dimensions
    dx0 = 0*x0
    dq0 = 0
    dz0 = 0*z0  #dummy variable analogous to xdot
    alg0 = 0*z0
    Nx = size(x0,1)
    Nq = size(q0, 1)
    Nz = size(z0,1)
    Nu = size(u0,1)
    #endregion

  # function to collect the variables in correct order for passing to DAE integrator
  function fun_DAE(Uout,DU,U,p,t)

    x = U[1:Nx]    #Doubt - how able to access Nx here. Need to make this general
    q = U[Nx+1]
    z = U[Nx+1+1 : Nx+1+Nz]
    u = p

    xdot = Plant_Eqns(x,z,u)[1]
    dq   = Plant_Eqns(x,z,u)[2]
    alg  = Plant_Eqns(x,z,u)[3]

  #Differential Equation for States (in form 0 = dx/dt - f(x,z,u))
    for nx = 1:Nx
      Uout[nx] = DU[nx] - xdot[nx]
    end

  #Differential Equation for Objective
    Uout[Nx+1] = DU[Nx+1] - dq[1]

  #Algebraic Equations  (in form 0 = g(x,z,u))
    for nz = 1:Nz
      Uout[Nx+1+nz] = alg[nz]
    end

  end

  ## Integrate the Plant for 1 time step (Solving DAE using Sundials)
  U₀ = vcat(x0,q0,z0)       #Initial Guess for all the variables (Differential, Objective, Algebraic) in DAE
  DU₀ = vcat(dx0,dq0,alg0)  #Initial Guess for their gradients (zero for the algebraic equations)
  differential_vars = vcat( ones(Nx+1), zeros(Nz))  #Identifier for which variables are differential and which are algebraic

  prob = DAEProblem(fun_DAE, DU₀,  U₀, tspan,  u0, differential_vars=differential_vars)

  #Solve the DAE
  sol = Sundials.solve(prob,IDA())
  xf = sol.u[end][1:Nx]
  qf = sol.u[end][Nx+1]
  zf = sol.u[end][Nx+2 : Nx+1+Nz]

    # # Plot the Profile for all [x, q, z]
    # plotly()
    # plot(sol)
  

  return xf, qf, zf
end

#test DAE Integrator
# Integrate_Plant_DAE(x0, q0, z0, u0, tspan)

##todo Create function for integrating ODE plant
# function Integrate_Plant_ODE(x0, q0, No_z0, u0, tspan)

# # end



##! Trying to solve for Steady State ->Not working
  #region

  # function fun_SS!(t, U)

  #   dU = 0*U
  #   Nx = size(x0,1)
  #   Nz = size(z0,1)
  #   # Nu = size(u0,1)
    
  #   x = U[1:2]    
  #   q = U[3:3]
  #   z = U[4:4]
  #   u0 = [0.1; 0.0]
  
  #   xdot = Plant_Eqns(x, z, u0)[1]
  #   dq   = Plant_Eqns(x, z, u0)[2]
  #   alg  = Plant_Eqns(x, z, u0)[3]


  #   for nx = 1:Nx
  #     dU[nx] = xdot[nx]
  #   end
  #   dU[Nx+1] = dq[1]
  #   for nz = 1:Nz
  #     dU[Nx+1+nz] = alg[nz]
  #   end

  #   return dU
  # end

  # U0 = vcat(x0, q0, z0)

  # fun_SS!(0.0, U0)

  # SS_prob = SteadyStateProblem(fun_SS!, U0)
  # SS_sol = solve(SS_prob, SSRootfind() )

  #endregion


