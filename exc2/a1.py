from time import sleep

import math
import numpy as np
import matplotlib.pyplot as plt

# define physical parameters
T0 = 20000
l0 = 10**-22
alpha = 10
beta = -0.5  
nH = 1
kb = 1.38 * 10**-16

# function lambda
def lmbd(T):
  if T <= T0:
    return l0 * (T/T0) ** alpha
  else: 
    return l0 * (T/T0) ** beta

# temporal differential of the temperature.
def function(T,t):
  return -2/(3 * kb) * nH * lmbd(T)

# runge kutta integrater of 2nd order. 
def rk(h,f,y,t):
  k1=f(y,T)
  k2=f(y+h*k1,t+h)
  return y + h/2*k1 + h/2*k2

# define a simulation class
class Simulation(object):
  def __init__(self, h, TInit, TStopp, maxErr):
    self.h=h
    self.t = 0.0 # time passed since the beginning.
    self.T=TInit
    self.TNew = 0.0
    self.THalve = 0.0
    self.err = 0.0 # estimated error of current step
    self.maxErr = maxErr # maximum allowed error per step
    self.i = 0 # number of time steps passed. 
    self.TStopp = TStopp
    
    # changing the size of the temperature array frequently is slowly.
    # therefore use a cache array. 
    self.cache = np.array(np.zeros(100)) # cash array.
    self.temp = np.array(()) # storage for the temperature developement.
    self.logTemp = np.array(()) # log of temperature.

    # same thing for the time
    self.tCache = np.array(np.zeros(100)) # cash for time array 
    self.timeArray = np.array(()) # storage for the time developement.
    
  def simulate(self):
    while self.T > self.TStopp:

      while True:
        # simulate with two h: one step with h and two steps with h/2
        self.TNew = rk(self.h,function,self.T,self.t) # simulation with normal stepsiz

        self.THalve = rk(self.h/2,function,self.T,self.t) # simulation with small stepsize 
        self.THalve = rk(self.h/2,function,self.THalve,self.t)

        self.err = abs(self.TNew - self.THalve) # get an estimate for the error in the current step.

        if self.err > self.maxErr:
          self.h = self.h/2
        elif 10 * self.err < self.maxErr:
          self.h = self.h*2
          break
        else:
          break

      self.T = self.THalve
      self.t += self.h

      # write simulation data in cache array.
      self.cache[self.i%100] = self.T
      self.tCache[self.i%100] = self.t

      # write the cache into the target array every 100 timesteps.
      if self.i%100 == 99:
        self.temp = np.append(self.temp, self.cache)
        self.timeArray = np.append(self.timeArray, self.tCache)

      self.i += 1


    self.temp = np.append(self.temp, self.cache[:self.i%100])
    self.timeArray = np.append(self.timeArray, self.tCache[:self.i%100])

    self.logTemp = np.zeros(len(self.temp))
    for j in range(len(self.logTemp)):
      if (self.temp[j] < 0.01):
        self.temp[j] = 0.01
      self.logTemp[j] = math.log(self.temp[j]) 

  def simulateFixedStep(self):
    while self.T > self.TStopp:
      self.t += self.h
      self.T = rk(self.h,function,self.T,self.t)

      # write simulation data in cache array.
      self.cache[self.i%100] = self.T
      self.tCache[self.i%100] = self.t

      # write the cache into the target array every 100 timesteps.
      if self.i%100 == 99:
        self.temp = np.append(self.temp, self.cache)
        self.timeArray = np.append(self.timeArray, self.tCache)

      self.i += 1


    self.temp = np.append(self.temp, self.cache[:self.i%100])
    self.timeArray = np.append(self.timeArray, self.tCache[:self.i%100])

    self.logTemp = np.zeros(len(self.temp))
    for j in range(len(self.logTemp)):
      if (self.temp[j] < 0.01):
        self.temp[j] = 0.01
      self.logTemp[j] = math.log(self.temp[j]) 

  def printResult(self):
    print "time passed: " + str(self.t)
    print "number of timesteps: " + str(self.i)

  def exportResult(self):
    return (self.logTemp, self.timeArray)

# parameters of the simulation
TStopp = 6000
h = 10. ** 10
TInit = 10.0 ** 7
T = TInit
t = 0.0
errMax = 10.0
errMaxSmall = 0.1

plt.subplot(111)

simVar = Simulation(h*10**4,TInit,TStopp,errMax)
simVar.simulate()
simVar.printResult()
(TVar,tVar) = simVar.exportResult()
plotVar = plt.plot(tVar,TVar, label = "Variable Stepsize\nerrMax = 10")

simVarSmall = Simulation(h,TInit,TStopp,errMaxSmall)
simVarSmall.simulate()
simVarSmall.printResult()
(TVarSmall,tVarSmall) = simVarSmall.exportResult()
plotVarSmall = plt.plot(tVarSmall,TVarSmall, label = "Variable Stepsize\nerrMax = 0.1")


simFixed = Simulation(h,TInit,TStopp,errMax)
simFixed.simulateFixedStep()
simFixed.printResult()
(TFixed,tFixed) = simFixed.exportResult()
plotFixed = plt.plot(tFixed,TFixed, label = "fixed Stepsize")

plt.legend(loc=1, borderaxespad=0.) #legend(bbox_to_anchor=(0, 1.05, 1, 0.1), loc=3, borderaxespad=0.)

plt.savefig("test.png")
plt.show((plotVar,plotVarSmall,plotFixed))

