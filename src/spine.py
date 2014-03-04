import sys, os,random
from walker import Walker 

class Spine:
	"""docstring for Spine"""
	def __init__(self,neck_width, D = 1, dt = 0.01):
		from math import sqrt
		# Must find a smart way to give diffusion "constant"
		# Must implement variable diffusion "constant"
		some_factor = 0.1
		self.neck_width = some_factor*neck_width
		self.neck_length = 1+random.random()*some_factor*self.neck_width
		self.head_height = 1-self.neck_length
		self.head_width = 0.5*(1-self.neck_width)
		self.spike_probability = 0.0003	# This number should probably be chosen more carefully
		self.Ions = []					# linked list of walkers
		self.dendrite_boundary = []		# linked list of walkers that have moved into the dendrite
		self.a = self.head_height/self.head_width
		self.drift = 0.0

		self._x1 = self._y1 = 1
		self._x0 = self._y0 = 0
		self.steplength = sqrt(2*D*dt)

	def SetDrift(self,val):
		self.drift = val

	def right_limit(self,x):
		return self.a*x-self.a*0.5*(1+self.neck_width)+self.neck_length

	def left_limit(self,x):
		return 1-self.a*x
	
	def AddDrift(self,Ion):
		Ion.r[1] -= self.drift

	def Solve(self):
		# Make some reasonable way to solve RW stuff
		if random.random()<self.spike_probability:
			self.AddSpike()
		for walker in self.Ions:
			walker.r[random.randint(0,1)] += (-1)**(random.randint(0,1))*self.steplength
			self.AddDrift(walker)
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
		if Ion.r[1]>self._y1:
			Ion.r[1] -= Ion.r[1]-self._y1
		elif Ion.r[0]>0.5 and Ion.r[1]<self.right_limit(Ion.r[0]):
			Ion.r[1] += self.right_limit(Ion.r[0])-Ion.r[1]
		elif Ion.r[0]<0.5 and Ion.r[1]<self.left_limit(Ion.r[0]):
			Ion.r[1] += self.left_limit(Ion.r[0])-Ion.r[1]
		elif Ion.r[1]<=self.neck_length:
			if Ion.r[0]<0.5*(1-self.neck_width):
				Ion.r[0] += 0.5*(1-self.neck_width)-Ion.r[0]
			elif Ion.r[0]>0.5*(1-self.neck_width):
				Ion.r[0] += 0.5*(1-self.neck_width)-Ion.r[0]
		if Ion.r[1]<self._y0:
			#ion has left spine
			self.dendrite_boundary.append(Ion)
			self.Ions.remove(Ion)
		
		# if Ion.r[0]>self._x1:
		# 	Ion.r[0] -= Ion.r[0]-self._x1
		# elif Ion.r[0]<self._x0:
		# 	Ion.r[0] -= Ion.r[0]-self._x0
		# if Ion.r[1]>self._y1:
		# 	Ion.r[1] -= Ion.r[1]-self._y1