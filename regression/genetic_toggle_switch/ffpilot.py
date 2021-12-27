#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

from lma.regression.models.genetic_toggle_switch.gtsRegression import FFPilotGTSRegression

class FFPilotGTSRegressionLMES(FFPilotGTSRegression):
    """subclass that allows for easy overriding of various parameters at point of testing
    """
    # def buildBasins(self):
    #     self.parser['both-basins'] = True
    #     return super(FFPilotGTSRegressionLMES, self).buildBasins()

    def _buildDefaultSimulationParameterDict(self):
        if self.parser['quick']:
            return {
                'errorGoal': .99,
                'errorGoalConfidence': .01,
                'ffpilotPilotOutput': False,
                'ffpilotPhaseOutput': False,
                'ffpilotStageOutputRaw': True,
                'ffpilotStageOutputSummary': True,
                'productionStageCountMinimum': 1e1,
                'binSpecies': True
                # 'condenseOutput': True,
                # 'writeInterval': 4.0,
                # 'pilotStageCount': 1e2,
                # 'writeInitialTrajectoryState': True,
                # 'writeFinalTrajectoryState': True,
                # 'writeLimitTracking': False,
            }
        else:
            # return {
            #     'errorGoal': .99,
            #     'errorGoalConfidence': .01,
            #     'ffpilotPilotOutput': False,
            #     'ffpilotPhaseOutput': False,
            #     'ffpilotStageOutputRaw': True,
            #     'ffpilotStageOutputSummary': True,
            #     'ffpilotMinimizeCost': False,
            #     'stepsPerWorkUnitPart': 4e5, #75000, #40000, #37000,
            #     'pilotStageCount': 1,
            #     'phaseZeroSamplingMultiplier': 5000,
            #     'productionStageCountMinimum': 1,
            #     'writeInitialTrajectoryState': False,
            #     'writeFinalTrajectoryState': False,
            #     'writeInterval': None,
            #     'writeLimitTracking': True,
            # }

            return {
                'errorGoal': .05,
                'errorGoalConfidence': .95,
                'ffpilotPilotOutput': True,
                'ffpilotPhaseOutput': False,
                'ffpilotStageOutputRaw': True,
                'ffpilotStageOutputSummary': True,
                'binSpecies': True,
                # 'condenseOutput': False,
                # 'pilotStageCount': 1e4,
                # 'writeInterval': 4.0,
                # 'writeInitialTrajectoryState': True,
                # 'writeFinalTrajectoryState': True,
                # 'writeLimitTracking': True,
            }

if __name__=='__main__':
    regression = FFPilotGTSRegressionLMES()
    regression.main()
