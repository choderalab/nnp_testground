
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
from openmm import app
import torch
from openmm import unit
torch._C._jit_set_nvfuser_enabled(False)

# Load the PDB file
pdb = app.PDBFile('ejm_54.pdb')

# Load the AMBER parameter and coordinate files
prmtop = app.AmberPrmtopFile('ejm_54.prm7')
inpcrd = app.AmberInpcrdFile('ejm_54.rst7')


ml_atoms = [idx for idx in range(4657, 4693)]
print(ml_atoms)
temp = 300 
dt = 0.0005  
step = 100_000 

potential = MLPotential('ani2x')
mm_system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
                             constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)


system = potential.createMixedSystem(prmtop.topology, mm_system, ml_atoms)

integrator = VerletIntegrator(dt * unit.picoseconds)
#platform = Platform.getPlatformByName("CUDA")
#prop = dict(CudaPrecision="mixed")
simulation = Simulation(prmtop.topology, system, integrator)#, platform, prop)
simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)

# # Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

simulation.reporters.append(DCDReporter(f'ani2x_simulation.dcd', 50, enforcePeriodicBox=True))
simulation.reporters.append(
    StateDataReporter(
        f'ani2x_simulation.csv',
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
