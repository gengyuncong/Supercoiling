#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

from lma.regression.models.self_regulating_gene.srgRegression import FFPilotSRGRegression

class FFPilotSRGRegressionLMES(FFPilotSRGRegression):
    """subclass that allows for easy overriding of various parameters at point of testing
    """
    def buildBasins(self):
        self.parser['both-basins'] = True
        return super(FFPilotSRGRegressionLMES, self).buildBasins()

    def buildTilings(self, **kwargs):
        # override passed value of `reverse`, if any
        kwargs['reverse'] = True

        return super(FFPilotSRGRegressionLMES, self).buildTilings(**kwargs)

    def _buildDefaultSimulationParameterDict(self):
        if self.parser['quick']:
            return {
                'ffpilotPilotOutput': True,
                'ffpilotPhaseOutput': True,
                'ffpilotStageOutputRaw': True,
                'ffpilotStageOutputSummary': True,
                'errorGoal': .9,
                'errorGoalConfidence': .01,
                'pilotStageCount': 10,
                'productionStageCountMinimum': 10,
                'writeInterval': 0.1,
                
                # 'batchSize': 1,
                # 'phaseZeroSamplingMultiplier': 1e2,
                # 'ffpilotMinimizeCost': True,
                # 'writeInitialTrajectoryState': False,
                # 'writeFinalTrajectoryState': False,
                # 'writeInterval': None,
                # 'writeLimitTracking': True,
                # 'stepsPerWorkUnitPart': 1e8,
            }
        else:
            return {
                'errorGoal': .2,

                'pilotStageCount': 1e4,
                'productionStageCountMinimum': 1e2,
                'writeInterval': 1.0,

                'ffpilotPilotOutput': True,
                'ffpilotPhaseOutput': False,
                'ffpilotStageOutputRaw': True,
                'ffpilotStageOutputSummary': True,

                # 'batchSize': 100,
                # 'errorGoalConfidence': .95,
                # 'ffpilotPilotOutput': True,
                # 'ffpilotPhaseOutput': False,
                # 'ffpilotStageOutputRaw': True,
                # 'ffpilotStageOutputSummary': True,
                # 'phaseZeroSamplingMultiplier': 1,
                # 'ffpilotMinimizeCost': True,
                # 'writeLimitTracking': False,
            }

if __name__=='__main__':
    regression = FFPilotSRGRegressionLMES()
    regression.main()
