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

# Exercise 4

The program above utilises a number of Monte Carlo systems to simulate a decay experiment. I employed a multilevel inheritance approach given the iterative nature of the experiment; all calculations lead into the subsequent ones, so the class system should reflect that. Should the program not have been relatively small, this approach may have become untenable with many parameters and the possibility of inheritance confusion. As it stands, this method is appropriate.

The inertial to lab frame velocity calculations were derived from Lorentz transformation equations; spherical to cartesian coordinates were calculated to make this more straightforward and to aid in the subsequent 3D and 2D density plots. Dimensional manipulation proved to be the most prevalent aspect of code. Due to the many different sorts of calculations, np.shape() was used extensively to test array sizes and then reshape arrays for purpose. This was especially important for the matrix calculations as the data had to be reformed and zipped correctly. 

The mpl_toolkits.mplot3d package was of great importance for both presenting data and testing the 3d distribution of direction vectors at the decay vertex, particularly for ensuring uniform spherical distribution in the inertial frame. Once transformed to the lab frame, the decay products can be seen to have an increased z velocity. At this point I went back and changed the distribution in the inertial frame to be hemispherical. This was because we are interested in +z velocities, as the tracking stations are at +z coords. In order to save computing power and get data for the most amount of events, this was the natural solution. It also superseded the necessity to remove -z velocity events from calculations involving the track parameters later on.

To solve the matrix equation for track parameter reconstruction, the script.linalg library was utilised. The lstsq() function computes a least squares solution to the matrix equation. It does this by computing a vector such that the norm |b-Mx| is minimised. The routine uses the ‘gelsd’ lapack-driver. 
 
The problem is solved in three steps:
 (1) Reduce the coefficient matrix M to bidiagonal form with
     Householder transformations, reducing the original problem
     into a "bidiagonal least squares problem" (BLS)
 (2) Solve the BLS using a divide and conquer approach.
 (3) Apply back all the Householder transformations to solve
     the original least squares problem.

 The divide and conquer algorithm makes very mild assumptions about floating point arithmetic.
