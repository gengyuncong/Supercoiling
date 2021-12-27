#!/usr/bin/env python3
import os

from lma.src.main import Sim

if __name__=='__main__':
    try:
        os.unlink('ffpilotSimulationInput.sfile')
    except OSError:
        pass

    sim = Sim(fpath='ffpilotSimulationInput.sfile', mode='w')

    tags = frozenset((('stage', 'Custom'),))
    ffpilotSimulationInput = sim.ffpilotSimulationInputs.getDatum(tags)

    # setting the stage
    ffpilotStage = ffpilotSimulationInput.ffpilot_stages.add()

    ffpilotStage.basin_id = 1
    ffpilotStage.tiling_id = 1

    # setting a phase
    ffpilotPhase = ffpilotStage.ffpilot_phases.add()

    ffpilotPhase.phase_id = 4
    ffpilotPhase.tile_id = 4
    ffpilotPhase.basin_id = 1
    ffpilotPhase.tiling_id = 1

    # setting a phase limit
    ffpilotPhaseLimit = ffpilotPhase.ffpilot_phase_limit

    ffpilotPhaseLimit.stop_condition = 1 # TRAJECTORY_COUNT
    ffpilotPhaseLimit.events_per_trajectory = 1
    ffpilotPhaseLimit.trajectories_per_phase = 4

    ffpilotPhaseLimit.uvalue = int(1e4)

    # setting a start point
    start_point = ffpilotPhase.start_points.add()
    start_point.species_coordinates = [0,1,2,3,4,5,6]
    start_point.count = 1
    start_point.times = [0.0]

    # ffpilot_phase_output = ffpilotSimulationInput.ffpilot_phase_output_list.ffpilot_phase_outputs.add()
    #
    # ffpilot_phase_output.tiling_id = 0
    # ffpilot_phase_output.basin_id = 0
    # ffpilot_phase_output.ffpilot_phase_index = 4
    #
    # end_point = ffpilot_phase_output.successful_trajectory_end_points.add()
    # end_point.species_coordinates = [0,1,2,3,4,5,6]
    # end_point.count = 1
    # end_point.times = [0.0]

    sim.ffpilotSimulationInputs.wtf(mode='w', excludedFields=('tiling', 'pilot_stage'))
