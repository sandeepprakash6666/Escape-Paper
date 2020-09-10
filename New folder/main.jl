
include("Plant.jl")
include("OCP.jl")


##* Declaring Simulation Parameters

T0_sim = 0
Tf_sim = 60
dt_sim = 1
N_sim = (Tf_sim - T0_sim) / dt_sim
N_sim = convert(Int, N_sim)
N_OCP = 30

yₛₚ = [1.5032 * ones(convert(Int, N_sim / 3));     0.9951 * ones(convert(Int, N_sim / 3));    0 * ones(convert(Int, N_sim / 3))]
                                                                            yₛₚ = append!(yₛₚ, 0 * ones(convert(Int, N_OCP)))

# Starting Conditions for Plant
uk₋₁ = [0.3]
xk₀ = [1.0; 1.0]
qk₀ = [0.0]
zk₀ = [0.0]

    #region ->Creating array for plotting
    Nx = size(xk₀, 1)
    Nz = size(zk₀, 1)
    Nu = size(uk₋₁, 1)

    x_plot = NaN * ones(Nx, N_sim + 1)
    q_plot = NaN * ones(1,  N_sim + 1)
    z_plot = NaN * ones(Nz, N_sim + 1)
    u_plot = NaN * ones(Nu, N_sim)

    # First Point for Plotting
    x_plot[:,1] = copy(xk₀)
    q_plot[:,1] = copy(qk₀)
    z_plot[:,1] = copy(zk₀)
    #endregion


##* Simulate Plant with NMPC

for k in 1:N_sim
    # k = 19
    global uk₋₁, xk₀, qk₀, zk₀

    yₛₚk = yₛₚ[k:k - 1 + N_OCP]

    # Solve OCP
    uk₀ = Solve_OCP(xk₀, uk₋₁, yₛₚk)

    # Simulate Plant
    xk₁, qk₁, zk₁ = Integrate_Plant_DAE(xk₀, qk₀, zk₀, uk₀, (0.0, dt_sim))

    #region -> Updating values for next iteration and Plotting
    xk₀   = copy(xk₁)
    qk₀   = copy(qk₁)
    zk₀   = copy(zk₁)
    uk₋₁  = copy(uk₀)

    ##Results for Plotting
    x_plot[:,k + 1]   = copy(xk₁)
    q_plot[1,k + 1]   = copy(qk₁)
    z_plot[:,k + 1]   = copy(zk₁)
    u_plot[:,k]     = copy(uk₀)
    #endregion

    
end

##* Plot Results

# plotly()
gr()

t_plot = collect(T0_sim:dt_sim:Tf_sim)

p11 = plot(t_plot, x_plot[1,:], label="x1")
p11 = plot!(t_plot, x_plot[2,:], label="x2")
p11 = plot!(t_plot, yₛₚ[1:N_sim + 1], label="y_sp",
                                                linetype=:steppost, linestyle=:dash)

p12 = plot(t_plot[1:end - 1], u_plot[1,:], label="u",
                                                linetype=:steppost)
fig1 = plot(p11, p12, layout=(2, 1))




