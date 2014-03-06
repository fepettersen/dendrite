import random, numpy as np
from spine import Spine
from diff import Diff as Diffusion

class Dendrite:
	"""docstring for Dendrite"""
	def __init__(self,dt,mesh_points,max_spine_size = 5):
		self.m = mesh_points
		dx = 1.0/(self.m+1)
		self.a = dt/(2*dx*dx)
		
		self.PDE = Diffusion(dt,dx,0)

		self.spines  = []
		self.indeces = []

		self.max_spine_contact_point = max_spine_size
		self.spine_spike_probability_threshold = 1e-3
		self.conversion_factor = 0.05
	
	def SetSpineSpikeProbabilityThreshold(self,value):
		self.spine_spike_probability_threshold = value

	# def SetInitialCondition(self,U0,D,v=0):
	# 	dx = 1.0/(self.m+1)
	# 	dt = self.a*(2*dx*dx)
	# 	drift = v*dt/(2*dx)
	# 	self.Up = U0
	# 	self.U = np.zeros(np.shape(U0))
	# 	self.PDE.Assemble(D,self.a,0,drift)
	# 	self.PDE.Precondition(self.PDE.A)

	def SetInitialCondition(self,U0,r_m,c_m,r_l):
		dx = 1.0/(self.m+1)
		dt = self.a*(2*dx*dx)
		self.Up = U0
		self.U = np.zeros(np.shape(U0))
		self.PDE.Assemble(dt,dx,r_m,c_m,r_l)
		self.PDE.Precondition(self.PDE.A)

	def AddSpine(self,drift=0.01):
		"""m is the total number of mesh-points on the dendrite"""
		m0 = 0
		n0 = self.m+1
		m0 = random.randint(2,self.m-2)
		while n0>=self.m:
			n0 = random.randint(m0+1,m0+self.max_spine_contact_point)

		self.spines.append(Spine(n0-m0))
		self.indeces.append([m0,n0])
		self.spines[-1].SetDrift(drift)
		print "spine added at index %d, %d"%(m0,n0)

	def Combine(self,spine,m,n):
		DX = 1.0/(n-m+1)
		for element in spine.dendrite_boundary:
			self.U[m+int(round(element.r[0]/DX))] += self.conversion_factor

		tmp = 0
		del spine.dendrite_boundary[:]
		for j in xrange(n-m):
			tmp += self.U[m+j]


	def Soma(self,x,t):
		pass

	def Solve(self):
		# self.PDE.Source = self.Soma
		self.U = self.PDE.Solve(self.Up)
		i = 0
		steps = 200
		for spine in self.spines:
			if random.random()<self.spine_spike_probability_threshold:
				spine.AddSpike()
			for j in xrange(steps):
				spine.Solve()
			# if len(spine.Ions)>0:
			# 	print "spine #",i," has %d \"ions\" "%len(spine.Ions)
			m0 = self.indeces[i][0]
			n0 = self.indeces[i][1]
			self.Combine(spine,m0,n0)
			i+=1
		self.Up = self.U.copy() 