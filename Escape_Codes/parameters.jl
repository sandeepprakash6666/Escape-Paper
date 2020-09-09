
module MyModule

## Unscaled Initial POints/ Guesses
x0_us = [60.0]
z0_us = [60.0; 60.0; 60.0]
u0_us = [0.0; 0.0]



ls_x = [0.0]
us_x = [100.0]

ls_z = [0.0; 0.0; 0.0]
us_z = [100.0; 100.0; 100.0]

ls_u = [0.0; 0.0]
us_u = [1.0; 10000]


Nx = size(x0_us, 1)
Nz = size(z0_us, 1)
Nu = size(u0_us, 1)


##Scaled Points

x0 = (x0_us - ls_x) ./ (us_x - ls_x)
z0 = (z0_us - ls_z) ./ (us_z - ls_z)
u0 = (u0_us - ls_u) ./ (us_u - ls_u)

export x0,      z0,     u0, 
       Nx,      Nz,     Nu, 
       ls_x,    ls_z,   ls_u, 
       us_x,    us_z,   us_u

end