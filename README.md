# UoB_CompPhys_yr3
Computational physics exercises for 3rd year - University of Bristol 

# Exercise 1

A remote overhead camera at a football stadium is suspended by three cables attached to the roof. Each cable is fed from a motorised 
drum so the camera can be moved around by changing the lengths of the cables. The camera has a mass of 50 kg and the attachment points 
are all 90 m from the centre of the pitch forming an equilateral triangle in a horizontal plane. The camera is moved in a horizontal 
plane at constant depth (7m) below the attachment points. You may ignore the mass of the cables.

You should use a 3D coordinate system  (𝑥,𝑦,𝑧) , in which  (0,0,0)  is above the centre of the pitch and the attachment points are at :

𝑝0=(−45√3,−45) ,  𝑝1=(0,90) ,  𝑝2=(45√3,−45)  metres.

# Exercise 2

In this problem we will find eigenvalues of the 1D Schrödinger equation using numerical methods.

The time-independent Schrödinger equation in 1D can be written :

𝐻𝜓 = 𝐸𝜓

Where the Hamiltonian  𝐻  is given by

𝐻 = −ℏ^2/2𝑚 𝑑^2/𝑑𝑥^2 + 𝑉

In order to find numerical solutions, we can divide the spatial dimension into  𝑁  discrete points,  𝑥𝑖 , and evaluate  𝜓  at each one. 
Given this, equation  1  becomes a matrix equation, with  𝜓  an  𝑁 -dimensional vector, and  𝐻  an  (𝑁×𝑁)  matrix. We can then find the 
eigenvalues and eigenfunctions of the equation using numerical methods.

# Exercise 3

The 1D diffusion equation is :

∂𝑢/∂𝑡 = 𝑘 ∂^2u/∂𝑥^2
 
You should discretize this equation onto  𝑁𝑥  space points, with separation  Δ𝑥=ℎ , and into timesteps  Δ𝑡=𝜏 . In the equations below, I use 
subscript  𝑖  as a space index, and superscript  𝑛  for time indices.

Having discretized the problem, you should use the implicit finite difference equation, as discussed in lectures :

𝑢𝑛+1𝑖−𝑢𝑛𝑖𝜏=𝑘𝑢𝑛+1𝑖+1−2𝑢𝑛+1𝑖+𝑢𝑛+1𝑖−1ℎ2
 
This can be written in matrix form  𝑢𝑛=𝑀𝑢𝑛+1  using :

𝑢𝑛𝑖=−𝛼𝑢𝑛+1𝑖−1+(1+2𝛼)𝑢𝑛+1𝑖−𝛼𝑢𝑛+1𝑖+1
 
where  𝛼=𝑘𝜏ℎ2 .

In the problems below, you are asked to solve the diffusion equation in the context of the heat equation. Here,  𝑘  is the thermal diffusivity, 
given by  𝑘=𝜆𝜌𝐶 , where  𝜆  is the thermal conductivity,  𝜌  is the density, and  𝐶  is the specific heat capacity. The questions below concern
an iron poker of length 50cm. You may take the thermal conductivity of iron to be a constant 59 W/m/K, its specific heat as 450 J/kg/K, and its 
density as 7,900 kg/m3. You can ignore heat loss along the length of the poker.
