import sys, os,random
from walker import Walker 

class Spine:
	"""docstring for Spine"""
	def __init__(self, D = 1, dt = 0.01):
		from math import sqrt
		# Must find a smart way to give diffusion "constant"
		# Must implement variable diffusion "constant"

		self.spike_probability = 0.0003	# This number should probably be chosen more carefully
		self.Ions = []					# linked list of walkers
		self.dendrite_boundary = []		# linked list of walkers that have moved into the dendrite

		self._x1 = self._y1 = 1
		self._x0 = self._y0 = 0
		self.steplength = sqrt(2*D*dt)

		
	def Solve(self):
		# Make some reasonable way to solve RW stuff
		if random.random()<self.spike_probability:
			self.AddSpike()
		for walker in self.Ions:
			walker.r[random.randint(0,1)] += (-1)**(random.randint(0,1))*self.steplength
			self.checkpos(walker)

	def AddSpike(self):
		# Release/spawn some number of ions at Post Synaptic Density
		amp = random.randint(1,15)		#should be drawn from Gaussian and have some meaning
		print "Spike of size ",amp
		for i in xrange(amp):
			#draw a random position in proximity of PSD
			self.Ions.append(Walker([random.random(),1-0.1*random.random()]))

	def checkpos(self,Ion):
		"""pos is a 2d vector with position [x,y]. 
		Implements reflecting boundaries and stores ions that have "left" the
		spine in the dendrite_boundary list."""
		if Ion.r[0]>self._x1:
			Ion.r[0] -= Ion.r[0]-self._x1
		elif Ion.r[0]<self._x0:
			Ion.r[0] -= Ion.r[0]-self._x0
		if Ion.r[1]>self._y1:
			Ion.r[1] -= Ion.r[1]-self._y1
		elif Ion.r[1]<self._y0:
			#ion has left spine
			self.dendrite_boundary.append(Ion)
			self.Ions.remove(Ion)