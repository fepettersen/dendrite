from dendrite import Dendrite
import numpy as np
import matplotlib.pyplot as mpl,matplotlib.animation as animation


n_spines = 4

m = 101
T = 100
dx = 1.0/(m+1)
dt = dx**2/2.0
D = np.ones(m)
x = np.linspace(0,1,m)

def InitialCondition(x,x0=0.1,sigma=0.07):
	return np.exp(-(x-x0)**2/(2*sigma**2))

dendrite = Dendrite(dt,m)
dendrite.SetInitialCondition(InitialCondition(x),D)

for i in xrange(n_spines):
	dendrite.AddSpine(drift=0.1)

im = []
fig = mpl.figure()
# dendrite.spines[0].AddSpike()
U_plot = np.zeros((m,2))
X,Y = np.meshgrid(np.linspace(0,1,m),np.linspace(0,dx,2))

for t in xrange(T):
	dendrite.Solve()
	im.append(mpl.plot(x,dendrite.U,'-b'))
	# U_plot[:,0] = dendrite.U.copy()
	# U_plot[:,1] = dendrite.U.copy()
	# im.append(mpl.contour(X,Y,U_plot))
	print "step %d of %d"%(t+1,T)

ani = animation.ArtistAnimation(fig,im,interval=180,blit=True)
mpl.show()