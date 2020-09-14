# using PlotlyJS
using Plots
include("OCP.jl")

dt = 1.0
NCP = 3

#region-> Create Models
##* Create Subproblem 1 

T01 = 0.0
Tf1 = 30.0
x0_us = [60.0]
Q_whb1 = vcat(1.0*ones(10,1), 1.2*ones(10,1), 0.8*ones(10,1)) *1.2539999996092727e6

t_plot1 = collect(T01:dt:Tf1) 

P1 = Build_OCP(x0_us, Q_whb1, Tf1 - T01)
# Naming variables to access solution
x1 = getindex(P1, :x)
T_tes1   = getindex(P1, :T_tes)
T_b1     = getindex(P1, :T_tes)
T_phb1   = getindex(P1, :T_phb)
T_whb1   = getindex(P1, :T_whb)
α1       = getindex(P1, :α)
Q_phb1   = getindex(P1, :Q_phb)
Constr_Coll_Diff0_1 = getindex(P1, :Constr_Coll_Diff0)
# z_ADMM1 = getindex(P1, :z_ADMM) 

JuMP.objective_sense(P1)
JuMP.objective_function(P1)


JuMP.set_NL_objective(P1, :Min, :(  $(Q_phb1[1])   ) )  #todo - Add if statements inside model to manually make subproblems



##* Create Subproblem 2
T02 = 30.0
Tf2 = 60.0
Q_whb2 = vcat(1.1*ones(10,1), 1.3*ones(10,1), 0.7*ones(10,1)) *1.2539999996092727e6 #todo - make into generic array or something

t_plot2 = collect(T02:dt:Tf2) 

P2 = Build_OCP(x0_us, Q_whb2, Tf2 - T02)
# Naming variables to access solution
x2 = getindex(P2, :x)
T_tes2   = getindex(P2, :T_tes)
T_b2     = getindex(P2, :T_tes)
T_phb2   = getindex(P2, :T_phb)
T_whb2  = getindex(P2, :T_whb)
α2       = getindex(P2, :α)
Q_phb2   = getindex(P2, :Q_phb)
Constr_Coll_Diff0_2 = getindex(P2, :Constr_Coll_Diff0)
z_ADMM2 = getindex(P2, :z_ADMM) 

#endregion

#region-> Solve Models
##* Solve SP1

JuMP.fix(z_ADMM1, 0.0; force = true)

optimize!(P1)
JuMP.termination_status(P1)
JuMP.solve_time(P1::Model)
JuMP.objective_value(P1)

    #extract solution to Julia variables
    star_T_tes1 = JuMP.value.(T_tes1[:,NCP])
                star_T_tes1 = cat(x0_us[1], star_T_tes1, dims = 1)     
    star_T_b1    = JuMP.value.(T_b1[:, NCP])
    star_T_phb1  = JuMP.value.(T_phb1[:, NCP])
    star_T_whb1  = JuMP.value.(T_whb1[:, NCP])

    star_α1      = JuMP.value.(α1[:])
    star_Q_phb1  = JuMP.value.(Q_phb1[:])

    star_z_ADMM1 = JuMP.value.(z_ADMM1)
    lambda1 = JuMP.dual.(Constr_Coll_Diff0_1)

##* Solve SP2 

JuMP.fix(z_ADMM2, 0.0; force = true)

optimize!(P2)
JuMP.termination_status(P2)
JuMP.solve_time(P2::Model)
JuMP.objective_value(P2)

    #extract solution to Julia variables
    star_T_tes2 = JuMP.value.(T_tes2[:,NCP])
                star_T_tes2 = cat(x0_us[1], star_T_tes2, dims = 1)     #todo - scaling to be made automatic
    star_T_b2    = JuMP.value.(T_b2[:, NCP])
    star_T_phb2  = JuMP.value.(T_phb2[:, NCP])
    star_T_whb2  = JuMP.value.(T_whb2[:, NCP])

    star_α2      = JuMP.value.(α2[:])
    star_Q_phb2  = JuMP.value.(Q_phb2[:])

    star_z_ADMM2 = JuMP.value.(z_ADMM2)
    lambda2 = JuMP.dual.(Constr_Coll_Diff0_2)


#endregion

##region-> Plotting
    #choose backend for plots
    plotly()
    # gr()

##* Plot Subproblem 1
    #Differential States
    p11 = plot(t_plot1, star_T_tes1)
    #Manipulated Variables
    p12 = plot(t_plot1[1:end-1], star_α1,
                                                                        marker = true, linetype = :steppost)
    p13 = plot(t_plot1[1:end-1], star_Q_phb1,        
                                                                        marker = true, linetype = :steppost)
    #Algebraic States
    p14 = plot(t_plot1[2:end], star_T_b1, 
                                                                        marker = true)
    p14 = plot!(t_plot1[2:end], star_T_phb1) 
                                                                        #marker = false)
    p15 = plot(t_plot1[2:end], star_T_whb1) 
                                                                # marker = false)

    # fig_P1 = plot(p11, p12, p13, p14, p15, layout = (5, 1))
    # gui(fig1)
    #endregion


##* Plot Subproblem 2
        #Differential States
        p21 = plot(t_plot2, star_T_tes2                  )
        #Manipulated Variables
        p22 = plot(t_plot2[1:end-1], star_α2,
                                                                            marker = true, linetype = :steppost)
        p23 = plot(t_plot2[1:end-1], star_Q_phb2,
                                                                            marker = true, linetype = :steppost)
        #Algebraic States
        p24 = plot(t_plot2[2:end], star_T_b2, 
                                                                            marker = true)
        p24 = plot!(t_plot2[2:end], star_T_phb2) 
                                                                            #marker = false)
        p25 = plot(t_plot2[2:end], star_T_whb2) 
                                                                    # marker = false)

    # fig_P2 = plot(p21, p22, p23, p24, p25, layout = (5, 1))


# fig_3 = plot(fig_P1, fig_P2, layout = (1,2))

fig_4 = plot(p11, p21)


##* Create Central Problem
T0 = T01
Tf = Tf2

Q_whb = vcat(Q_whb1, Q_whb2)

t_plot = collect(T0:dt:Tf) 

Centr = Build_OCP(x0_us, Q_whb, Tf - T0)
# Naming variables to access solution
x = getindex(Centr, :x)
T_tes   = getindex(Centr, :T_tes)
T_b     = getindex(Centr, :T_tes)
T_phb   = getindex(Centr, :T_phb)
T_whb   = getindex(Centr, :T_whb)
α       = getindex(Centr, :α)
Q_phb   = getindex(Centr, :Q_phb)
Constr_Coll_Diff0 = getindex(Centr, :Constr_Coll_Diff0)
z_ADMM = getindex(Centr, :z_ADMM) 


##* Solve Central Problem

JuMP.fix(z_ADMM, 0.0; force = true)

optimize!(Centr)
JuMP.termination_status(Centr)
JuMP.solve_time(Centr::Model)
JuMP.objective_value(Centr)

#extract solution to Julia variables
star_T_tes = JuMP.value.(T_tes[:,NCP])
             star_T_tes = cat(x0_us[1], star_T_tes, dims = 1)     #todo - scaling to be made automatic
star_T_b    = JuMP.value.(T_b[:, NCP])
star_T_phb  = JuMP.value.(T_phb[:, NCP])
star_T_whb  = JuMP.value.(T_whb[:, NCP])

star_α      = JuMP.value.(α[:])
star_Q_phb  = JuMP.value.(Q_phb[:])

star_z_ADMM = JuMP.value.(z_ADMM)
lambda = JuMP.dual.(Constr_Coll_Diff0)

##* 