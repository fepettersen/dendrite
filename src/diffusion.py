import numpy as np, time
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.animation as animation
from math import sqrt

def Assemble_and_Decompose(D, a, b):
	#Assemble the matrix A which will be constant as long as dt is constant, 
	#and make a LU decomposition of A
	m,n = np.shape(D)
	N = m*n
	A = np.zeros((N,N))
	k=0
	for i in xrange(m):
		A[i,i+n] = -2*a*(D[0,i]+D[1,i]);
		A[N-1-i,N-i-n-1] = -2*a*(D[m-1,m-1-i]+D[m-2,m-1-i]);
		for j in xrange(n):
			if j==0:
				A[k,k+1] = -2*b*(D[i,j+1]+D[i,j])
				if i==0:
					A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
				elif i==m-1:
					A[k,k-1] = 0
					A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
				else:
					A[k,k-1] = 0
					A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
			elif j==n-1:
				A[k,k-1] = -2*b*(D[i,j-1]+D[i,j])
				if i==0:
					A[k,k+1] = 0
					A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j-1]
				elif i==m-1:
					A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) +2*b*D[i,j-1]
				else:
					A[k,k+1] = 0
					A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) +2*b*D[i,j-1]
			elif j!=0:
				A[k,k+1] = -b*(D[i,j+1]+D[i,j])
				A[k,k-1] = -b*(D[i,j-1]+D[i,j])
				if i==0:
					A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
				elif i==m-1:
					A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
				else:
					A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
			if (k>(m-1) and k<(N-m)):
				A[k,k+m] = -a*(D[i+1,j]+D[i,j]);
				A[k,k-m] = -a*(D[i-1,j]+D[i,j]);
			k+=1
	return A;
def Assemble1d(D,a):
	n = len(D)
	N = n**2
	A = np.zeros((n,n))
	k=0
	A[k,k] = 1+2*a*(D[1]+D[0])
	A[k,k+1] = -2*a*(D[0]+D[1])
	k+=1
	for i in xrange(1,n-1):
		A[k,k] = 1+a*D[i+1]+2*a*D[i] +a*D[i-1]
		A[k,k+1] = -a*(D[i+1]+D[i])
		A[k,k-1] = -a*(D[i-1]+D[i])
		k+=1
	A[k,k] = 1+2*a*(D[n-2]+D[n-1])
	A[k,k-1] = -2*a*(D[n-1]+D[n-2])
	return A

def Precondition(A):
	n = int(sqrt(np.shape(A)[0]))
	N = n**2
	H = []
	D = []
	Aa = []
	i=0
	a = np.zeros((n,n))
	b = np.zeros((n,n))
	c = np.zeros((n,n))

	b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
	c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

	D.append(np.linalg.inv(b))
	H.append(-1*np.dot(D[-1],c))
	for i in xrange(1,n-1):
		a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

		Aa.append(a.copy())
		D.append(np.linalg.inv(b+np.dot(a,H[-1])))
		H.append(-1*np.dot(D[-1],c))
		
	i=n-1

	a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
	Aa.append(a.copy())
	b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
	D.append(np.linalg.inv(b+np.dot(a,H[-1])))

	return H,D,Aa
	
def modified_blockTridiag(H,D,a,Up):
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


def polyreg(x,y,m):
	deg=m
	A = np.zeros((deg,deg))
	b = np.zeros(deg)
	for i in xrange(deg):
		tempe = 0
		for l in xrange(m):
			b[i] += x[l]**i*y[l]
			tempe += x[l]**i*y[l]
		for j in xrange(i,deg):
			tmp = 0
			for k in xrange(m):
				tmp += x[k]**(i+j)
			A[i,j] = A[j,i] = tmp
	a = np.linalg.solve(A,b)

	for l in xrange(m):
		y[l] =0
		for i in xrange(deg):
			y[l] += a[i]*x[l]**i
	return y

def Gaussian(x,y,t,x0 = 0.5,y0 = 1.0, sigma_x = 0.65, sigma_y = 0.05):
	return np.cos(np.pi*t)**2*np.exp(-(x-x0)**2/(2*sigma_x**2)-(y-y0)**2/(2*sigma_y**2))


def f(x,y,t=0):
	return np.exp(-np.pi*np.pi*t)*np.cos(np.pi*x)*np.cos(np.pi*y)
def S(x,y,t):
	pi = np.pi
	py =pi*y
	px = pi*x
	# return pi*np.exp(-t*pi*pi)*(-pi*np.cos(px)*np.cos(py)*(-2*x-2*y+pi)+\
	# 	np.sin(px)*np.cos(py) +np.sin(py)*np.cos(px));
	return pi*np.exp(-t*pi*pi)*(2*pi*np.cos(px)*np.cos(py)*(x+y-0.5) \
		+ np.cos(px)*np.sin(py) +np.sin(px)*np.cos(py));

d = 2
M = 26
Hc = 100
T = 200
dt = 0.001
x,y = np.meshgrid(np.linspace(0,1,M),np.linspace(0,1,M))
# D = x+y
D = np.ones((M,M))*0.5
dx = 1./(M-1)
a = dt/(2*dx*dx)
b = dt/(2*dx*dx)
# A = Assemble_and_Decompose(D,a,b)
if d==1:
	# A = Assemble1d(np.linspace(0,1,M),a)
	x = np.linspace(0,1,M)
	A = Assemble1d(np.ones(M),a)
	x0 = 3; x1 = 10
	Up = f(x,np.zeros(M))
	im = []
	fig = mpl.figure()
	for t in xrange(T):
		im.append(mpl.plot(np.linspace(0,1,M),Up,'b-'))
		U = np.linalg.solve(A,Up)
		U[x0:x1] += (1/np.sqrt(Hc))*np.random.random_sample(x1-x0)
		U[:] = polyreg(x,U,M)
		Up = U.copy()
	ani = animation.ArtistAnimation(fig,im)
	mpl.show()
else:
	A = Assemble_and_Decompose(D,a,b)
	H,D,a = Precondition(A)
	# np.savetxt('delete_me.txt',A)
	# for i in xrange(M**2):
	# 	for j in xrange(M**2):
	# 		print A[i,j]," ",
	# 	print ""
	# for i in xrange(M**2):
	# 	print A[i,i]
	# print "a = ",a
	Up = 0.4*Gaussian(x,y,0,sigma_x=0.5,sigma_y=0.3,x0 = 0.5,y0 = 0.75)
	mpl.ion()
	fig = mpl.figure()
	ax = fig.add_subplot(111,projection="3d")


	for t in xrange(T):
		wframe = ax.plot_wireframe(x,y,Up)
		ax.set_xlabel('x axis')
		ax.set_ylabel('y axis')
		if t==0:
			ax.set_autoscaley_on(False)
			# pass
		mpl.draw()
		ax.collections.remove(wframe)
		Up += dt*Gaussian(x,y,t*dt)
		Uptmp = Up.reshape(-1)
		Utmp = modified_blockTridiag(H,D,a,Uptmp)
		# Utmp = blockTridiag(A,Uptmp)
		U = Utmp.reshape(M,M)
		Up = U.copy()

	# fig2 = mpl.figure()
	# ax2 = fig2.add_subplot(111,projection="3d")
	# blergh = ax2.plot_wireframe(x,y,D)
	# raw_input('press enter>>')