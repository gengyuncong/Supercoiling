#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

from lma.regression.models.genetic_toggle_switch.gtsRegression import ReplicateGTSRegression

class ReplicateGTSRegressionLMES(ReplicateGTSRegression):
    """subclass that allows for easy overriding of various parameters at point of testing
    """
    def _buildDefaultSimulationParameterDict(self):
        if self.parser['quick']:
            return {'maxTime': 1e1,
                    'stepsPerWorkUnitPart': 1e8,
                    # 'binSpecies': True,
                    # 'condenseOutput': False
                    'writeInterval': 1e0,
                    'degreeAdvancementWriteInterval': 1e0,
                    'orderParameterWriteInterval': 1e0,
                    }
        else:
            return {'maxTime': 1e5 + 3.0,
                    'stepsPerWorkUnitPart': 1e8,
                    # 'binSpecies': True,
                    # 'condenseOutput': False,
                    'writeInterval': 4.0,
                    # 'degreeAdvancementWriteInterval': 4.0,
                    # 'orderParameterWriteInterval': 1e1,
                    'writeInitialTrajectoryState': True,
                    'writeFinalTrajectoryState': True,
                    }

if __name__=='__main__':
    regression = ReplicateGTSRegressionLMES()
    regression.main()