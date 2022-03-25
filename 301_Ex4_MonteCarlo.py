# beam simulation functions

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import numpy.random as rand
import scipy.linalg as linalg

rand.seed(42)

NUM_PARTICLES = 1000 # particles/decay events
MU_VELOCITY = 2000 # m s^-1 . Mean injection velocity in lab frame
SIGMA_VELOCITY = 50 # velocity standard deviation
TAU_MEAN = 2.5e-3 # mean lifetime of particle
DAUGHTER_SPEED = 5.0e+3 # m s^-1 . isotropic decay speed in inertial frame
STATIONS = (30, 35, 40, 45) # m . Station distances
DETECTOR_RESOLUTION = 0.01 # m . Hit resolution ry,rx.

class BeamGeneration:

    def __init__(self, NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN):
        self.N = NUM_PARTICLES
        self.mu = MU_VELOCITY
        self.sigma = SIGMA_VELOCITY
        self.tau = TAU_MEAN

    def beam_velocity(self):

        zeros = np.zeros((self.N, 3))

        beam = rand.normal(self.mu, self.sigma, self.N)
        beam = np.reshape(beam,(self.N,1))

        # Insert 1D array of speeds into beam velocity vector arrays
        indices = np.full((self.N,1),2)
        np.put_along_axis(zeros,indices,beam,1)
        return zeros

    def decay_time(self):

        # Exponential decay of particles, returns distribution array.
        t = rand.exponential(self.tau, self.N)
        return t

    def decay_vertex(self):
        
        # Simple d = v*t to return array of decay vertices.
        vertex = self.beam_velocity() * np.reshape(self.decay_time(),(self.N,1))
        return vertex

    def plot_generation(self):
        
        fig, axs = plt.subplots(1, 3, figsize=(12,5))
        fig.suptitle("Beam Generation Test")
        axs[0].title.set_text('Beam Velocity')
        axs[1].title.set_text('Decay Time')
        axs[2].title.set_text('Decay Vertex')

        n_bins = 50

        bv = self.beam_velocity()
        dt = self.decay_time()
        dv = self.decay_vertex()
        
        # Slice only relevant columns of arrays for histogram plot.
        axs[0].hist(bv[:,2], bins=n_bins)
        axs[1].hist(dt, bins=n_bins)
        axs[2].hist(dv[:,2], bins=n_bins)

# decay simulation functions

class BeamDecay(BeamGeneration):

    def __init__(self, NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED):
        super().__init__(NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN)
        self.d_speed = DAUGHTER_SPEED
        
    def decay_direction_inertial(self):
        
        # Generate dierction angles, only +ve z as tracking stations are +ve z. 
        phi = rand.uniform(0, 2*np.pi, self.N)
        costheta = rand.uniform(0, 1, self.N)
        theta = np.arccos(costheta)

        # Transform to cartesian coordinates, makes Lorentz velocity addition easier.
        x_hat = np.sin(theta) * np.cos(phi)
        y_hat = np.sin(theta) * np.sin(phi)
        z_hat = np.cos(theta)

        cartesian = np.dstack((x_hat,y_hat,z_hat))

        return cartesian
    
    def decay_direction_lab(self):
        
        inert_dir = self.decay_direction_inertial()
        v = self.beam_velocity()
        v = np.expand_dims(v, axis=0)
        c = 3.0e+8

        # Lorentz factors in z-direction.
        beta = v[0,:,2] / c
        gamma = 1 / np.sqrt(1 - beta**2)
        gamma = np.expand_dims(gamma, axis=0)
        gamma = np.expand_dims(gamma, axis=2)

        # Normalise velocities to get direction vectors.
        inert_norm = np.linalg.norm(inert_dir, axis=2)
        inert_norm = np.expand_dims(inert_norm, axis=2)
        r_hat = inert_dir / inert_norm
        u = r_hat * self.d_speed

        # Calculate dot product of frame velocity and particle velocity.
        u_dot_v = (v*u).sum(axis=2)
        u_dot_v = np.expand_dims(u_dot_v, axis=2)

        # Calculate Lorentz velocity vector addition
        lab_vel = (1/(1+(u_dot_v/(c**2)))) * (u/gamma + v + (gamma/((c**2)*(gamma+1)))*u_dot_v*v)

        # Again, normalise.
        lab_norm = np.linalg.norm(lab_vel, axis=2)
        lab_norm = np.expand_dims(lab_norm, axis=2)
        lab_dir = lab_vel / lab_norm

        return lab_dir
    
    def plot_inertial(self):

        fig = plt.figure()
        fig.suptitle("Inertial Frame Decay Direction")
        ax = Axes3D(fig)

        inert_dir = self.decay_direction_inertial()
        ix = inert_dir[0,:,0]
        iy = inert_dir[0,:,1]
        iz = inert_dir[0,:,2]

        ax.scatter(iz, ix, iy, color='red', s=1)
        ax.set_xlim(-1,1)
        plt.show()

    def plot_lab(self):

        fig = plt.figure()
        fig.suptitle("Lab Frame Decay Direction")
        ax = Axes3D(fig)

        lab_dir = self.decay_direction_lab()
        lx = lab_dir[0,:,0]
        ly = lab_dir[0,:,1]
        lz = lab_dir[0,:,2]

        ax.scatter(lz, lx, ly, color='red', s=1)
        ax.set_xlim(-1,1)
        plt.show()

# tracking functions

class DaughterPropagation(BeamDecay):

    def __init__(self, NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED, STATIONS):
        super().__init__(NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED)
        self.stations = STATIONS

    def track_parameters(self):

        vertex = self.decay_vertex()
        direction = self.decay_direction_lab()

        # Adjust injection distance so that all decays occur before x-y intercepts.
        inject_to_z0 = np.full((self.N),50)
        decay_to_z0 = inject_to_z0 - vertex[:,2]

        # Calculate dx/dz and dy/dz gradients
        grad_x = direction[0,:,0] / direction[0,:,2]
        grad_y = direction[0,:,1] / direction[0,:,2]

        # Calculate x-y coords of decay product at z=0
        intercept_x = decay_to_z0 * grad_x
        intercept_y = decay_to_z0 * grad_y
        
        parameters = grad_x, grad_y, intercept_x, intercept_y
        return parameters

    def station_hits(self):

        parameters = self.track_parameters()

        gradient = np.array([parameters[0], parameters[1]])
        intercept = np.array([parameters[2], parameters[3]])
        
        # Calculate hit station coordinates.
        hit = [[],[],[],[]]
        for i in range(len(self.stations)):
            hit[i] = intercept + (self.stations[i] * gradient)

        hit_coords = np.transpose(hit,(0,2,1))
        return hit_coords

    def plot_hits(self):

        hit_data = self.station_hits()

        max = math.ceil(np.amax(abs(hit_data[3,:,:])))
        fig = plt.figure()
        fig.suptitle("Tracking Station Hits")
        
        ax = Axes3D(fig, rect=(-0.5,-0.5,2,2))
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, 0.4, 0.4, 1]))
        
        ax.scatter(hit_data[0,:,0], hit_data[0,:,1], zs=self.stations[0], zdir='x', s=1, c='red')
        ax.scatter(hit_data[1,:,0], hit_data[1,:,1], zs=self.stations[1], zdir='x', s=1, c='orangered')
        ax.scatter(hit_data[2,:,0], hit_data[2,:,1], zs=self.stations[2], zdir='x', s=1, c='orange')
        ax.scatter(hit_data[3,:,0], hit_data[3,:,1], zs=self.stations[3], zdir='x', s=1, c='gold')
                
        station_delta = (self.stations[1] - self.stations[0])/2
        ax.set_xlim(self.stations[0]-station_delta, self.stations[3]+station_delta)
        ax.set_ylim(-max, +max)
        ax.set_zlim(-max, +max)
        
        plt.show()

# smearing function

class HitSmearing(DaughterPropagation):

    def __init__(self, NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED, STATIONS, DETECTOR_RESOLUTION):
        super().__init__(NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED, STATIONS)
        self.resolution = DETECTOR_RESOLUTION

    def position_smear(self):

        hit_res = self.resolution

        smear_x = [[],[],[],[]]
        smear_y = [[],[],[],[]]
        for i in range(len(self.stations)):
            smear_x[i] = rand.normal(0, hit_res, self.N)
            smear_y[i] = rand.normal(0, hit_res, self.N)

        smear = [smear_x, smear_y]
        smear = np.transpose(smear,(1,2,0))

        return smear
        
    def measured_hits(self):

        true_hits = self.station_hits()
        smear = self.position_smear()
        #print(np.shape(true_hits))
        #print(np.shape(smear))

        hit = true_hits + smear
        return hit

    def smear_distribution_plot(self):

        true_hits = self.station_hits()
        smeared_hits = self.measured_hits()

        max = math.ceil(np.amax(abs(true_hits[3,:,:])))

        fig, axs = plt.subplots(1, 2, figsize=(11,5))
        fig.suptitle("Hit Comparison (Station 1)")
        axs[0].title.set_text('True')
        axs[1].title.set_text('Measured')
        
        axs[0].set_ylim(-max, +max)
        axs[1].set_xlim(-max, +max)
        
        axs[0].hist2d(true_hits[0,:,0], true_hits[0,:,1], bins=70)
        axs[1].hist2d(smeared_hits[0,:,0], smeared_hits[0,:,1], bins=70)

# track reconstruction functions

class TrackReconstruction(HitSmearing):

    def __init__(self, NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED, STATIONS, DETECTOR_RESOLUTION):
        super().__init__(NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED, STATIONS, DETECTOR_RESOLUTION)

    def matrix(self):
        
        sub_matrix = [[],[],[],[]]
        for i in range(len(self.stations)):
            sub_matrix[i] = np.array([
                [self.stations[i], 0, 1, 0],
                [0, self.stations[i], 0, 1]])

        matrix = np.reshape(sub_matrix,(-1, 4))
        return matrix

    def solve_matrix(self):

        matrix = self.matrix()
        hits = self.measured_hits()
        hits = np.transpose(hits, (1,0,2))
        hits = np.reshape(hits, (self.N,-1), 'A')
        hits = np.transpose(hits)
       
        cs = linalg.lstsq(matrix, hits)
        
        return cs

# vertex reconstruction functions

class VertexReconstruction(TrackReconstruction):

    def __init__(self, NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED, STATIONS, DETECTOR_RESOLUTION):
        super().__init__(NUM_PARTICLES, MU_VELOCITY, SIGMA_VELOCITY, TAU_MEAN, DAUGHTER_SPEED, STATIONS, DETECTOR_RESOLUTION)

    def reconstruct_vertex(self):

        cs = self.solve_matrix()
        parameters = cs[0]

        p = np.transpose(parameters)

        z = -(p[:,0]*p[:,2] + p[:,1]*p[:,3]) / (p[:,0]**2 + p[:,1]**2)

        inj_to_dec = z + 50
        
        return inj_to_dec

    def vertex_plot(self):
        
        true_vertex = self.decay_vertex()
        reconstructed_vertex = self.reconstruct_vertex()

        print("Reconstructed Decay Vertex Mean: "+str(np.mean(reconstructed_vertex)))
        print("True Decay Vertex Mean: "+str(np.mean(true_vertex[:,2])))

        fig, axs = plt.subplots(1, 2, figsize=(11,5))
        fig.suptitle("Decay Vertex Comparison")
        axs[0].title.set_text('True')
        axs[1].title.set_text('Reconstructed')
        
        axs[0].hist(true_vertex[:,2], bins=50)
        axs[1].hist(reconstructed_vertex, bins=50) 