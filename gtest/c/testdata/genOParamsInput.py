#!/usr/bin/env python3
import os

from lma.src.main import Sim

if __name__=='__main__':
    try:
        os.unlink('oparamsInput.sfile')
    except OSError:
        pass
    sim = Sim(fpath='oparamsInput.sfile', mode='w')

    tags = frozenset([('prefix', 'Test')])
    simulationInput = sim.simulationInputs.getDatum(tags)

    # setting the oparams
    oparams = simulationInput.order_parameters.order_parameters

    oparamLinear = oparams.add()
    oparamLinear.id = 0
    oparamLinear.type = 0

    oparamLinear.species_ids = [0,1,2,3,4,5]
    oparamLinear.species_coefficients = [-1, -2, -2, 1, 2, 2]

    oparamExponential = oparams.add()
    oparamExponential.id = 1
    oparamExponential.type = 100

    oparamExponential.species_ids = [0,1,2,3,4,5]
    oparamExponential.species_coefficients = [-1, -2, -2, 1, 2, 2]
    oparamExponential.species_exponents = [2,1,1,2,1,1]

    oparamBasins = oparams.add()
    oparamBasins.id = 2
    oparamBasins.type = 1000

    for basinCounts in ((4,16,1,0,0,0,0), (0,0,0,4,16,1,0)):
        basin = oparamBasins.basins.add()
        basin.species_count = basinCounts

    sim.simulationInputs.wtf(mode='w')  #, excludedFields=('tiling', 'pilot_stage'))
