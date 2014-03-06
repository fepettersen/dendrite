from dendrite import Dendrite
from spine import Spine
import numpy as np, glob
import matplotlib.pyplot as mpl,matplotlib.animation as animation


n_spines = 0

m = 101
T = 100
dx = 1.0/(m+1)
dt = dx**2/2.0
D = np.ones(m)
drift = 0.3

Rm = 10**4
Cm = 0.002
radius = 1e-6
Rl = 0.352


rm = Rm/(2*np.pi*radius)
cm = Cm*2*np.pi*radius
rl = Rl/(np.pi*radius**2)

x = np.linspace(0,1,m)

def InitialCondition(x,x0=0.1,sigma=0.07):
	return np.exp(-(x-x0)**2/(2*sigma**2))
	# return np.zeros(np.shape(x))

dendrite = Dendrite(dt,m)
# dendrite.SetInitialCondition(InitialCondition(x),D*2,drift)
dendrite.SetInitialCondition(InitialCondition(x),rm,cm,rl)

for i in xrange(n_spines):
	dendrite.AddSpine(drift=0.1)
dat_spine = Spine(4)
dat_spine.AddSpike()

"""
for i in xrange(50):
	dat_spine.Solve()
	outfile = open('RW_verify_%03d.txt'%i,'w')
	outfile.write('%d\n'%len(dat_spine.Ions))
	for j in xrange(len(dat_spine.Ions)):
		pos = dat_spine.Ions[j].r
		outfile.write('%g %g \n'%(pos[0],pos[1]))
	outfile.close()
print "a = ",dat_spine.a
mpl.ion()
for step in sorted(glob.glob('RW_verify_*.txt')):
	tmp = open(step,'r')
	nwalkers = int(tmp.readline())
	for i in xrange(nwalkers):
		pos = tmp.readline().split()
		mpl.plot(float(pos[0]),float(pos[1]),'r-x')
		mpl.draw()

mpl.plot(x,dat_spine.left_limit(x),'r-')
mpl.hold('on')
mpl.plot(x,dat_spine.right_limit(x),'b-')
mpl.plot(x,np.ones(np.shape(x))*dat_spine.neck_length)
mpl.show()
"""
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
	# print "step %d of %d"%(t+1,T)

ani = animation.ArtistAnimation(fig,im,interval=180,blit=True)
mpl.show()
