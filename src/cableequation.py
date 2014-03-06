import numpy as np

class CableEquation:
	"""docstring for CableEquation"""
	def __init__(self,DT,dx,dy,x0=0,x1=1):
		self.preconditioned = False
		self.t = 0
		self.dt = DT
		self.m = int(1.0/dx)-1
		self.n = n = 1
		self.x = np.linspace(x0,x1,self.m)

	def Assemble(self,dt,dx,rm,cm,rl):
		"""Assemble the matrix A which is used to solve the linear system at each time-step"""
		tau = cm*rm
		lambda2 = rm/rl
		# print "tau = ",tau,"\nlambda^2 = ",lambda2
		m = self.m
		n = 1
		N = m*n
		self.A = np.zeros((N,N))
		if n==1:
			k=0
			self.A[k,k] = 1 + dt/tau + (2*lambda2*dt)/(tau*dx*dx)
			self.A[k,k+1] = -2*(lambda2*dt)/(tau*dx*dx)
			k+=1
			for i in xrange(1,m-1):
				self.A[k,k] = 1 + dt/tau + (2*lambda2*dt)/(tau*dx*dx)
				self.A[k,k+1] = -(lambda2*dt)/(tau*dx*dx)
				self.A[k,k-1] = -(lambda2*dt)/(tau*dx*dx)
				k+=1
			self.A[k,k] = 1 + dt/tau + (2*lambda2*dt)/(tau*dx*dx)
			self.A[k,k-1] = -2*(lambda2*dt)/(tau*dx*dx)

	def Precondition(self,A):
		self.preconditioned = True
		n = self.n
		m = self.m
		# n = int(np.sqrt(np.shape(A)[0]))
		N = m*n
		self.H = []
		self.D_ = []
		self.Aa = []
		i=0
		a = np.zeros((n,n))
		b = np.zeros((n,n))
		c = np.zeros((n,n))

		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

		self.D_.append(np.linalg.inv(b))
		self.H.append(-1*np.dot(self.D_[-1],c))
		for i in xrange(1,m-1):
			a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
			b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
			c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

			self.Aa.append(a.copy())
			self.D_.append(np.linalg.inv(b+np.dot(a,self.H[-1])))
			self.H.append(-1*np.dot(self.D_[-1],c))
			
		i=m-1

		a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
		self.Aa.append(a.copy())
		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		self.D_.append(np.linalg.inv(b+np.dot(a,self.H[-1])))

	def modified_Block_Tridiag(self,H,D,a,Up):
		n = self.n
		m = self.m
		N = m*n
		g = []
		k = np.zeros(n)
		x = np.zeros(N)


		i = 0
		k[:] = Up[i*n:(i+1)*n]
		g.append(np.dot(D[0],k))
		for i in xrange(1,m-1):
			gtmp = g[-1]
			k[:] = Up[i*n:(i+1)*n]
			g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))
		i = m-1
		# print "Up[%d:%d] = "%(i*n,(i+1)*n),Up[i*n],",",Up[(i+1)*n]
		gtmp = g[-1]
		k[:] = Up[i*n:(i+1)*n]
		g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))

		x[i*n:(i+1)*n] = g[i]
		for i in xrange(m-2,-1,-1):
			x[i*n:(i+1)*n] = g[i]+np.dot(H[i],x[(i+1)*n:(i+2)*n])

		return x

	def Source(self,x,t):
		return 0

	def Solve(self,Up):
		Up += self.dt*self.Source(self.x,self.t)
		if self.preconditioned:
			U = self.modified_Block_Tridiag(self.H,self.D_,self.Aa,Up)
		else:
			U = self.Block_Tridiag(Up)
		self.t += self.dt
		return U