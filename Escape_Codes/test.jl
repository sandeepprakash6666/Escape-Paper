

using Revise
includet("MyFile.jl")
using Main.Mymodule


# include("MyFile.jl")
# using Main.Mymodule

#accesssing from workspace directly
println(x)  
println(y)  #Does not update it's value
mygreeting()


println(Main.Mymodule.x)
println(Main.Mymodule.y)    #Does not update it's value
println(Main.Mymodule.mygreeting())

# includet("C:/Users/sandeepp/NTNU/o365_Sandeep's_teams - General/PhD sem 1/Escape Paper/Escape_Codes/parameters.jl")



