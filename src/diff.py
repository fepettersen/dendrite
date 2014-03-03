import numpy as np

class Diff:
	"""docstring for Diff"""
	def __init__(self,DT):
		self.preconditioned = False
		self.t = 0
		self.dt = DT

	def Assemble(D,a,b):
		#Assemble the matrix A which will be constant as long as dt is constant, 
		#and make a LU decomposition of A
		if not np.shape(D):
			import sys
			print "Cannot assemble from scalar diffusion constant. Use D = np.ones(..)*D. \n exiting \n "
			sys.exit(1)

		m,n = np.shape(D)
		N = m*n
		self.A = np.zeros((N,N))
		k=0
		for i in xrange(m):
			self.A[i,i+n] = -2*a*(D[0,i]+D[1,i]);
			self.A[N-1-i,N-i-n-1] = -2*a*(D[m-1,m-1-i]+D[m-2,m-1-i]);
			for j in xrange(n):
				if j==0:
					self.A[k,k+1] = -2*b*(D[i,j+1]+D[i,j])
					if i==0:
						self.A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
					elif i==m-1:
						self.A[k,k-1] = 0
						self.A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
					else:
						self.A[k,k-1] = 0
						self.A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
				elif j==n-1:
					self.A[k,k-1] = -2*b*(D[i,j-1]+D[i,j])
					if i==0:
						self.A[k,k+1] = 0
						self.A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j-1]
					elif i==m-1:
						self.A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) +2*b*D[i,j-1]
					else:
						self.A[k,k+1] = 0
						self.A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) +2*b*D[i,j-1]
				elif j!=0:
					self.A[k,k+1] = -b*(D[i,j+1]+D[i,j])
					self.A[k,k-1] = -b*(D[i,j-1]+D[i,j])
					if i==0:
						self.A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
					elif i==m-1:
						self.A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
					else:
						self.A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
				if (k>(m-1) and k<(N-m)):
					self.A[k,k+m] = -a*(D[i+1,j]+D[i,j]);
					self.A[k,k-m] = -a*(D[i-1,j]+D[i,j]);
				k+=1

	def Precondition(self,A):
		self.preconditioned =
		n = int(sqrt(np.shape(A)[0]))
		N = n**2
		H = []
		D = []
		self.Aa = []
		i=0
		a = np.zeros((n,n))
		b = np.zeros((n,n))
		c = np.zeros((n,n))

		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

		self.D_.append(np.linalg.inv(b))
		self.H.append(-1*np.dot(self.D_[-1],c))
		for i in xrange(1,n-1):
			a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
			b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
			c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

			self.Aa.append(a.copy())
			self.D_.append(np.linalg.inv(b+np.dot(a,self.H[-1])))
			self.H.append(-1*np.dot(self.D_[-1],c))
			
		i=n-1

		a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
		self.Aa.append(a.copy())
		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		self.D_.append(np.linalg.inv(b+np.dot(a,self.H[-1])))

	def modified_Block_Tridiag(self,H,D,a,Up):
		n = int(sqrt(np.shape(A)[0]))
		N = n**2
		g = []
		k = np.zeros(n)
		x = np.zeros(N)

		i = 0
		k[:] = Up[i*n:(i+1)*n]
		g.append(np.dot(D[0],k))
		for i in xrange(1,n-1):
			gtmp = g[-1]
			k[:] = Up[i*n:(i+1)*n]
			g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))
		i = n-1
		gtmp = g[-1]
		k[:] = Up[i*n:(i+1)*n]
		g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))

		x[i*n:(i+1)*n] = g[i]
		for i in xrange(n-2,-1,-1):
			x[i*n:(i+1)*n] = g[i]+np.dot(H[i],x[(i+1)*n:(i+2)*n])

		return x

	def Source(self,x,t):
		pass

	def Solve(self,Up):
		Up += dt*self.Source(x,self.t):
		if precondition:
			U = self.modified_Block_Tridiag(self.H,self.D_,self.Aa,Up)
		else:
			U = self.Block_Tridiag(Up)
		self.t += self.dt
		return U