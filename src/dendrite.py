
class Dendrite:
	"""docstring for Dendrite"""
	def __init__(self):
		PDE = Diffusion()
	
	def SetInitialCondition(self,U0):
		self.Up = U0
		self.U = np.zeros(np.shape(U0))
		PDE.Assemble(U0)
		PDE.Precondition(PDE.A)

	def AddSpine(self,x0,x1):
		m0,n0 = self.MapAreaToIndex(x0,x1)
		self.spines.append(Spine())

	def MapAreaToIndex(self,x0,x1):
		m = int(round(x0/dx))
		n = int(round(x1/dx))
		return m,n

	def Combine(self,tmp,m,n):
		pass

	def Solve(self):
		self.U = PDE.solve(self.Up)
		for spine in self.spines:
			tmp = spine.solve()
			self.Combine(self.U,tmp,m,n)