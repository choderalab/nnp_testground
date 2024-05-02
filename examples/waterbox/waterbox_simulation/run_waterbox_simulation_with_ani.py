
from openmm.app import (
    StateDataReporter,
    DCDReporter,
    CharmmPsfFile,
    Simulation,
    CharmmCrdFile,
)
from openmm import (
    unit,
    VerletIntegrator,
    Platform,
)
from openmmml import MLPotential
import torch

torch._C._jit_set_nvfuser_enabled(False)

# read in topology and coordinates
psf=CharmmPsfFile("tip125.psf")
crd=CharmmCrdFile("tip125_cptequil.crd")

# system specific parameters
temp = 300 
dt = 0.0005  
step = 100_000 

# set up box, potential and simulation object
box_size = 15.5554 #smaller box size for 125 water molecules
psf.setBox(box_size * unit.angstrom, box_size * unit.angstrom, box_size * unit.angstrom)
potential = MLPotential('ani2x')
system = potential.createSystem(psf.topology)
integrator = VerletIntegrator(dt * unit.picoseconds)
platform = Platform.getPlatformByName("CUDA")
prop = dict(CudaPrecision="mixed")
simulation = Simulation(psf.topology, system, integrator, platform, prop)
simulation.context.setPositions(crd.positions)
simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)

# # Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

simulation.reporters.append(DCDReporter(f'ani2x_waterbox_simulation.dcd', 50, enforcePeriodicBox=True))
simulation.reporters.append(
    StateDataReporter(
        f'ani2x_waterbox_simulation.csv',
        reportInterval=50,
        step=True,
        time=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        density=True,
        speed=True,
        separator="\t",
    )
)

if step > 0:
    print("\nMD run: %s steps" % step)
    simulation.step(step)
