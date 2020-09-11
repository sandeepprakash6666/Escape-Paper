

include("OCP.jl")

dt = 1.0
NCP = 3

##* Create Subproblem 1 
#region
x0_us = [60.0]
Q_whb1 = vcat(1.2*ones(10,1), ones(10,1), 0.8*ones(10,1)) *1.2539999996092727e6

T01 = 0.0
Tf1 = 30.0
t_plot1 = collect(T01:dt:Tf1) 

P1 = Build_OCP(x0_us, Q_whb1, Tf1)

optimize!(P1)
JuMP.termination_status(P1)
JuMP.solve_time(P1::Model)
JuMP.objective_value(P1)


    T_tes1   = getindex(P1, :T_tes)
    T_b1     = getindex(P1, :T_tes)
    T_phb1   = getindex(P1, :T_phb)
    T_whb1   = getindex(P1, :T_whb)
    α1       = getindex(P1, :α)
    Q_phb1   = getindex(P1, :Q_phb)


    star_T_tes1 = JuMP.value.(T_tes1[:,NCP])
                star_T_tes1 = cat(x0[1]*100, star_T_tes1, dims = 1)     #todo - scaling to be made automatic
    star_T_b1    = JuMP.value.(T_b1[:, NCP])
    star_T_phb1  = JuMP.value.(T_phb1[:, NCP])
    star_T_whb1  = JuMP.value.(T_whb1[:, NCP])

    star_α1      = JuMP.value.(α1[:])
    star_Q_phb1  = JuMP.value.(Q_phb1[:])


##  Plotting        
t_plot1      
    #choose backend for plots
    # plotly()
    gr()

        #Differential States
        p11 = plot(t_plot1, star_T_tes1,                label = "T_tes1")
        #Manipulated Variables
        p12 = plot(t_plot1[1:end-1], star_α1,           label = "α1",
                                                                            marker = true, linetype = :steppost)
        p13 = plot(t_plot1[1:end-1], star_Q_phb1,        label = "Q_phb1 kJ/hr",
                                                                            marker = true, linetype = :steppost)
        #Algebraic States
        p14 = plot(t_plot1[2:end], star_T_b1,             label = "T_b1", 
                                                                            marker = true)
        p14 = plot!(t_plot1[2:end], star_T_phb1,          label = "T_phb1") 
                                                                            #marker = false)
        p15 = plot(t_plot1[2:end], star_T_whb1,          label = "T_whb1") 
                                                                    # marker = false)

    fig1 = plot(p11, p12, p13, p14, p15, layout = (5, 1))
#endregion

##* Create Subproblem 2

Q_whb2 = vcat(1.0*ones(10,1), 1.3*ones(10,1), 0.7*ones(10,1)) *1.2539999996092727e6


T02 = 30.0
Tf2 = 60.0
t_plot2 = collect(T02:dt:Tf2) 

P2 = Build_OCP(x0_us, Q_whb2, Tf1)

optimize!(P2)
JuMP.termination_status(P2)
JuMP.solve_time(P2::Model)
JuMP.objective_value(P2)


    T_tes2   = getindex(P2, :T_tes)
    T_b2     = getindex(P2, :T_tes)
    T_phb2   = getindex(P2, :T_phb)
    T_whb2  = getindex(P2, :T_whb)
    α2       = getindex(P2, :α)
    Q_phb2   = getindex(P2, :Q_phb)


    star_T_tes2 = JuMP.value.(T_tes2[:,NCP])
                star_T_tes2 = cat(x0[1]*100, star_T_tes2, dims = 1)     #todo - scaling to be made automatic
    star_T_b2    = JuMP.value.(T_b2[:, NCP])
    star_T_phb2  = JuMP.value.(T_phb2[:, NCP])
    star_T_whb2  = JuMP.value.(T_whb2[:, NCP])

    star_α2      = JuMP.value.(α2[:])
    star_Q_phb2  = JuMP.value.(Q_phb2[:])


                        
    #choose backend for plots
    # plotly()
    gr()

        #Differential States
        p11_2 = plot(t_plot2, star_T_tes2,                  label = "T_tes2")
        #Manipulated Variables
        p12_2 = plot(t_plot2[1:end-1], star_α2,             label = "α2",
                                                                            marker = true, linetype = :steppost)
        p13_2 = plot(t_plot2[1:end-1], star_Q_phb2,        label = "Q_phb2 kJ/hr",
                                                                            marker = true, linetype = :steppost)
        #Algebraic States
        p14_2 = plot(t_plot2[2:end], star_T_b2,             label = "T_b2", 
                                                                            marker = true)
        p14_2 = plot!(t_plot2[2:end], star_T_phb2,          label = "T_phb2") 
                                                                            #marker = false)
        p15_2 = plot(t_plot2[2:end], star_T_whb2,          label = "T_whb2") 
                                                                    # marker = false)

    fig2 = plot(p11_2, p12_2, p13_2, p14_2, p15_2, layout = (5, 1))

## s

plot(fig1, fig2, layout = (1,2))