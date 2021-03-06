# UoB_CompPhys_yr3
Computational physics exercises for 3rd year - University of Bristol 

# Exercise 1

A remote overhead camera at a football stadium is suspended by three cables attached to the roof. Each cable is fed from a motorised 
drum so the camera can be moved around by changing the lengths of the cables. The camera has a mass of 50 kg and the attachment points 
are all 90 m from the centre of the pitch forming an equilateral triangle in a horizontal plane. The camera is moved in a horizontal 
plane at constant depth (7m) below the attachment points. You may ignore the mass of the cables.

You should use a 3D coordinate system  (π₯,π¦,π§) , in which  (0,0,0)  is above the centre of the pitch and the attachment points are at :

π0=(β45β3,β45) ,  π1=(0,90) ,  π2=(45β3,β45)  metres.

# Exercise 2

In this problem we will find eigenvalues of the 1D SchrΓΆdinger equation using numerical methods.

The time-independent SchrΓΆdinger equation in 1D can be written :

π»π = πΈπ

Where the Hamiltonian  π»  is given by

π» = ββ^2/2π π^2/ππ₯^2 + π

In order to find numerical solutions, we can divide the spatial dimension into  π  discrete points,  π₯π , and evaluate  π  at each one. 
Given this, equation  1  becomes a matrix equation, with  π  an  π -dimensional vector, and  π»  an  (πΓπ)  matrix. We can then find the 
eigenvalues and eigenfunctions of the equation using numerical methods.

# Exercise 3

The 1D diffusion equation is :

βπ’/βπ‘ = π β^2u/βπ₯^2
 
You should discretize this equation onto  ππ₯  space points, with separation  Ξπ₯=β , and into timesteps  Ξπ‘=π . In the equations below, I use 
subscript  π  as a space index, and superscript  π  for time indices.

Having discretized the problem, you should use the implicit finite difference equation, as discussed in lectures :

π’π+1πβπ’πππ=ππ’π+1π+1β2π’π+1π+π’π+1πβ1β2
 
This can be written in matrix form  π’π=ππ’π+1  using :

π’ππ=βπΌπ’π+1πβ1+(1+2πΌ)π’π+1πβπΌπ’π+1π+1
 
where  πΌ=ππβ2 .

In the problems below, you are asked to solve the diffusion equation in the context of the heat equation. Here,  π  is the thermal diffusivity, 
given by  π=πππΆ , where  π  is the thermal conductivity,  π  is the density, and  πΆ  is the specific heat capacity. The questions below concern
an iron poker of length 50cm. You may take the thermal conductivity of iron to be a constant 59 W/m/K, its specific heat as 450 J/kg/K, and its 
density as 7,900 kg/m3. You can ignore heat loss along the length of the poker.

# Exercise 4

The program above utilises a number of Monte Carlo systems to simulate a decay experiment. I employed a multilevel inheritance approach given the iterative nature of the experiment; all calculations lead into the subsequent ones, so the class system should reflect that. Should the program not have been relatively small, this approach may have become untenable with many parameters and the possibility of inheritance confusion. As it stands, this method is appropriate.

The inertial to lab frame velocity calculations were derived from Lorentz transformation equations; spherical to cartesian coordinates were calculated to make this more straightforward and to aid in the subsequent 3D and 2D density plots. Dimensional manipulation proved to be the most prevalent aspect of code. Due to the many different sorts of calculations, np.shape() was used extensively to test array sizes and then reshape arrays for purpose. This was especially important for the matrix calculations as the data had to be reformed and zipped correctly. 

The mpl_toolkits.mplot3d package was of great importance for both presenting data and testing the 3d distribution of direction vectors at the decay vertex, particularly for ensuring uniform spherical distribution in the inertial frame. Once transformed to the lab frame, the decay products can be seen to have an increased z velocity. At this point I went back and changed the distribution in the inertial frame to be hemispherical. This was because we are interested in +z velocities, as the tracking stations are at +z coords. In order to save computing power and get data for the most amount of events, this was the natural solution. It also superseded the necessity to remove -z velocity events from calculations involving the track parameters later on.

To solve the matrix equation for track parameter reconstruction, the script.linalg library was utilised. The lstsq() function computes a least squares solution to the matrix equation. It does this by computing a vector such that the norm |b-Mx| is minimised. The routine uses the βgelsdβ lapack-driver. 
 
The problem is solved in three steps:
 (1) Reduce the coefficient matrix M to bidiagonal form with
     Householder transformations, reducing the original problem
     into a "bidiagonal least squares problem" (BLS)
 (2) Solve the BLS using a divide and conquer approach.
 (3) Apply back all the Householder transformations to solve
     the original least squares problem.

 The divide and conquer algorithm makes very mild assumptions about floating point arithmetic.
