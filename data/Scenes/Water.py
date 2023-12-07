import splishsplash as sph
import numpy as np

def init(base):
    global counter
    print("init test")
    counter = 1
    
def step():
    global counter
    sim = sph.Simulation.getCurrent()
    # air fluid model
    fluid = sim.getFluidModel(0) if sim.getFluidModel.getId() == "Air" else sim.getFluidModel(1)
    tm = sph.TimeManager.getCurrent()
    print(fluid.getPosition(0))
    print(tm.getTime())
    print(counter)
    counter += 1
    print("---")