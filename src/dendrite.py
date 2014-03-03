
class Dendrite:
	"""docstring for Dendrite"""
	def __init__(self,mesh_points,max_spine_size = 5):
		self.m = mesh_points
		self.max_spine_contact_point = max_spine_size
		PDE = Diffusion()
	
	def SetInitialCondition(self,U0):
		self.Up = U0
		self.U = np.zeros(np.shape(U0))
		PDE.Assemble(U0)
		PDE.Precondition(PDE.A)

	def AddSpine(self):
		"""m is the total number of mesh-points on the dendrite"""
		m0 = n0 = 0
		m0 = random.random_integer(2,self.m-2)
		while n0<m0 and (n0-m0)>self.max_spine_contact_point:
			n0 = random.random_integer(0,self.m)
		self.spines.append(Spine())
		self.indeces.append([m0,n0])

	def Combine(self,tmp,m,n):
		tmp = self.U[m:n]


	def Solve(self):
		PDE.Source = self.Soma
		self.U = PDE.solve(self.Up)
		i = 0
		for spine in self.spines:
			if random.random()<self.spine_spike_probability_threshold:
				spine.AddSpike()
			tmp = spine.solve()
			m0 = self.indeces[i][0]
			n0 = self.indeces[i][1]
			self.Combine(self.U,tmp,m0,n0)