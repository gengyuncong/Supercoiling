#!/usr/bin/env python3

from lma import Sim

if __name__=='__main__':
    sim = Sim(fpath='ffpilotSupervisorInput.sfile', mode='w')

    tags = frozenset((('stage', 'Custom'),))
    ffpilotSimulationInput = sim.ffpilotSimulationInputs[tags]

    # setting the stage
    ffpilotStage = ffpilotSimulationInput.ffpilot_stage_list.ffpilot_stages.add()

    ffpilotStage.basin_index = 0
    ffpilotStage.tiling_id = 0

    # setting a phase
    ffpilotPhase = ffpilotStage.ffpilot_phases.add()

    ffpilotPhase.ffpilot_phase_index = 1
    ffpilotPhase.tile_index = 1
    ffpilotPhase.basin_index = 0
    ffpilotPhase.tiling_id = 0

    # setting a phase limit
    ffpilotPhaseLimit = ffpilotPhase.ffpilot_phase_limit

    ffpilotPhaseLimit.stop_condition = 1 # TRAJECTORY_COUNT
    ffpilotPhaseLimit.events_per_trajectory = 1
    ffpilotPhaseLimit.trajectories_per_phase = 4

    ffpilotPhaseLimit.uvalue = int(1e4)

    # setting a start point
    start_point = ffpilotPhase.start_points.add()
    start_point.species_coordinates = [4, 8, 1, 0, 0, 0, 0]
    start_point.count = 1
    start_point.times = [0.0]

    sim.ffpilotSimulationInputs.wtf(mode='w')