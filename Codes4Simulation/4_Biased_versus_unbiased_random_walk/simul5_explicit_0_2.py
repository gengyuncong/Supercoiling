import sys
import math

###########################################
def wirte_head(f):
	f.write('<?xml version="1.0" encoding="UTF-8"?> \n')
	f.write('<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1"> \n')
	f.write('	<model id="This_is_a_model_for_FtsZ" substanceUnits="item" timeUnits="second" volumeUnits="litre" extentUnits="item"> \n')
	f.write('		<listOfUnitDefinitions> \n')
	f.write('			<unitDefinition id="per_second"> \n')
	f.write('				<listOfUnits> \n')
	f.write('					<unit kind="second" exponent="-1" scale="0" multiplier="1"/> \n')
	f.write('				</listOfUnits> \n')
	f.write('			</unitDefinition> \n')
	f.write('			<unitDefinition id="per_item"> \n')
	f.write('				<listOfUnits> \n')
	f.write('					<unit kind="item" exponent="-1" scale="0" multiplier="1"/> \n')
	f.write('				</listOfUnits> \n')
	f.write('			</unitDefinition> \n')
	f.write('			<unitDefinition id="per_item_per_second"> \n')
	f.write('				<listOfUnits> \n')
	f.write('					<unit kind="item"   exponent="-1" scale="0" multiplier="1"/> \n')
	f.write('					<unit kind="second" exponent="-1" scale="0" multiplier="1"/> \n')
	f.write('				</listOfUnits> \n')
	f.write('			</unitDefinition> \n')
	f.write('			<unitDefinition id="per_molar_per_second"> \n')
	f.write('				<listOfUnits> \n')
	f.write('					<unit kind="litre"  exponent="1" scale="0" multiplier="1"/> \n')
	f.write('					<unit kind="mole"   exponent="-1" scale="0" multiplier="1"/> \n')
	f.write('					<unit kind="second" exponent="-1" scale="0" multiplier="1"/> \n')
	f.write('				</listOfUnits> \n')
	f.write('			</unitDefinition> \n')
	f.write('		</listOfUnitDefinitions> \n')
	f.write('		<listOfCompartments> \n')
	f.write('			<compartment id="cell" size="8e-15" spatialDimensions="3" constant="true"/> \n')
	f.write('		</listOfCompartments> \n')
	f.write('		<listOfSpecies> \n') 


def write_sp(f, name, amount):
	f.write('			<species id="{}"'.format(name)+' compartment="cell" initialAmount="{}" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>\n'.format(amount))


def write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	if (not Products=="NA"):
		f.write('				<listOfProducts> \n')
		for i in Products:
			f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
		f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw> \n')
	f.write('					<math xmlns="http://www.w3.org/1998/Math/MathML"> \n')
	f.write('						<apply> \n')
	f.write('							<times/> \n')
	f.write('							<ci> rate </ci> \n')
	if (not Propensity_sp=="NA"):
		for i in Propensity_sp:
			f.write('							<ci> {} </ci> \n'.format(i))
	f.write('						</apply> \n')
	f.write('					</math> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="rate" value="{}" units="per_second"/> \n'.format(Propensity_para))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')

def write_reaction2(f, index, Reactants,Reactants_stoi, Products, Products_stoi, Propensity_sp, Propensity_para):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in range(0,len(Reactants)):
		f.write('					<speciesReference species="{}"'.format(Reactants[i]))
		f.write(' stoichiometry="{}" constant="true"/> \n'.format(Reactants_stoi[i]))
	f.write('				</listOfReactants> \n')
	if (not Products=="NA"):
		f.write('				<listOfProducts> \n')
		for i in range(0,len(Products)):
			f.write('					<speciesReference species="{}"'.format(Products[i]))
			f.write(' stoichiometry="{}" constant="true"/> \n'.format(Products_stoi[i]))
		f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw> \n')
	f.write('					<math xmlns="http://www.w3.org/1998/Math/MathML"> \n')
	f.write('						<apply> \n')
	f.write('							<times/> \n')
	f.write('							<ci> rate </ci> \n')
	if (not Propensity_sp=="NA"):
		for i in Propensity_sp:
			f.write('							<ci> {} </ci> \n'.format(i))
	f.write('						</apply> \n')
	f.write('					</math> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="rate" value="{}" units="per_second"/> \n'.format(Propensity_para))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction3(f, index, Reactants,Reactants_stoi, Products, Products_stoi, Propensity_sp, Propensity_para):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in range(0,len(Reactants)):
		f.write('					<speciesReference species="{}"'.format(Reactants[i]))
		f.write(' stoichiometry="{}" constant="true"/> \n'.format(Reactants_stoi[i]))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in range(0,len(Products)):
		f.write('					<speciesReference species="{}"'.format(Products[i]))
		f.write(' stoichiometry="{}" constant="true"/> \n'.format(Products_stoi[i]))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw> \n')
	f.write('					<math xmlns="http://www.w3.org/1998/Math/MathML"> \n')
	f.write('						<apply> \n')
	f.write('							<times/> \n')
	f.write('							<cn> {} </cn> \n'.format(Propensity_para))
	f.write('							<ci> {} </ci> \n'.format(Reactants[0]))
	f.write('							<apply> \n')
	f.write('								<power/> \n')
	f.write('								<ci> {} </ci> \n'.format(Reactants[1]))
	f.write('								<cn> {} </cn> \n'.format(Reactants_stoi[1]))
	f.write('							</apply> \n')
	f.write('						</apply> \n')
	f.write('					</math> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction_stall(f, index, Reactants, Products, lk0, Inf, PCoil_threshold, NCoil_threshold):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="RNAPStallPropensity"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="lk0" value="{}"/> \n'.format(lk0))
	f.write('						<localParameter id="torqueDomainBoundary0" value="0.0126"/> \n')
	f.write('						<localParameter id="torqueDomainBoundary1" value="0.0375"/> \n')
	f.write('						<localParameter id="torqueValue0" value="516.33"/> \n')
	f.write('						<localParameter id="torqueValue1" value="6.51"/> \n')
	f.write('						<localParameter id="torqueValue2" value="173.6"/> \n')
	f.write('                       <localParameter id="stallCutoff" value="10.5"/> \n')
#	f.write('						<localParameter id="stallCutoff" value="18.5"/> \n')
	f.write('						<localParameter id="stallRate" value="{}"/> \n'.format(Inf))
	f.write('						<localParameter id="ncoilCutoff" value="{}"/> \n'.format(NCoil_threshold))
	f.write('						<localParameter id="pcoilCutoff" value="{}"/> \n'.format(PCoil_threshold))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction_resume(f, index, Reactants, Products, lk0, resumeRate, PCoil_threshold, NCoil_threshold):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="RNAPResumePropensity"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="lk0" value="{}"/> \n'.format(lk0))
	f.write('						<localParameter id="torqueDomainBoundary0" value="0.0126"/> \n')
	f.write('						<localParameter id="torqueDomainBoundary1" value="0.0375"/> \n')
	f.write('						<localParameter id="torqueValue0" value="516.33"/> \n')
	f.write('						<localParameter id="torqueValue1" value="6.51"/> \n')
	f.write('						<localParameter id="torqueValue2" value="173.6"/> \n')
	f.write('                       <localParameter id="stallCutoff" value="10.5"/> \n')
#	f.write('						<localParameter id="stallCutoff" value="18.5"/> \n')
	f.write('						<localParameter id="resumeRate" value="{}"/> \n'.format(resumeRate))
	f.write('						<localParameter id="pcoilCutoff" value="{}"/> \n'.format(PCoil_threshold))
	f.write('						<localParameter id="ncoilCutoff" value="{}"/> \n'.format(NCoil_threshold))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction_initiation(f, index, Reactants, Products, Boundary0, Boundary1, initiationRate, lk0, basalLevel):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="InitiationPropensity"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="Boundary0" value="{}"/> \n'.format(Boundary0))
	f.write('						<localParameter id="Boundary1" value="{}"/> \n'.format(Boundary1))
	f.write('						<localParameter id="initiationRate" value="{}"/> \n'.format(initiationRate))
	f.write('						<localParameter id="lk0" value="{}"/> \n'.format(lk0))
	f.write('						<localParameter id="basalLevel" value="{}"/> \n'.format(basalLevel))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction_gyrase_catalysis(f, index, Reactants, Products, sigma_t, epsilon, catalyticRate, lk0):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="GyraseCatalysisPropensity"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="sigma_t" value="{}"/> \n'.format(sigma_t))
	f.write('						<localParameter id="epsilon" value="{}"/> \n'.format(epsilon))
	f.write('						<localParameter id="catalyticRate" value="{}"/> \n'.format(catalyticRate))
	f.write('						<localParameter id="lk0" value="{}"/> \n'.format(lk0))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction_TopoI_catalysis(f, index, Reactants, Products, Products_stoi, catalyticRate, lk0):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in range(0,len(Products)):
		f.write('					<speciesReference species="{}"'.format(Products[i]))
		f.write(' stoichiometry="{}" constant="true"/> \n'.format(Products_stoi[i]))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="TopoICatalysisPropensity"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="catalyticRate" value="{}"/> \n'.format(catalyticRate))
	f.write('						<localParameter id="lk0" value="{}"/> \n'.format(lk0))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction_SC_diffusion(f, index, Reactants, Products, Products_stoi, diffusion_rate):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in range(0,len(Products)):
		f.write('					<speciesReference species="{}"'.format(Products[i]))
		f.write(' stoichiometry="{}" constant="true"/> \n'.format(Products_stoi[i]))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="SCDiffusionPropensityNew"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="diffusionCoefficient" value="{}"/> \n'.format(diffusion_rate))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')

def write_reaction_stall_diffusion(f, index, Reactants, Products, diffusion_rate):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="SCDiffusionPropensity"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="diffusionCoefficient" value="{}"/> \n'.format(diffusion_rate))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')


def write_reaction_RNAP_translocation(f, index, Reactants, Products, Inf):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="RNAPtranslocation"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="translocationRate" value="{}"/> \n'.format(Inf))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')

def write_reaction_mRNA_degradation_elongation(f, index, Reactants, Products, RNase_elongation_rate):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="mRNA_degradation_elongation"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="RNase_elongation_rate" value="{}"/> \n'.format(RNase_elongation_rate))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')

def write_PI_reaction(f, index, Reactants, Products, Inf):
	f.write('			<reaction id="Reaction{}" reversible="false" fast="false"> \n'.format(index))
	f.write('				<listOfReactants> \n')
	for i in Reactants:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfReactants> \n')
	f.write('				<listOfProducts> \n')
	for i in Products:
		f.write('					<speciesReference species="{}" stoichiometry="1" constant="true"/> \n'.format(i))
	f.write('				</listOfProducts> \n')
	f.write('				<kineticLaw metaid="PIPropensity"> \n')
	f.write('					<listOfLocalParameters> \n')
	f.write('						<localParameter id="inhibitionRate" value="{}"/> \n'.format(Inf))
	f.write('					</listOfLocalParameters> \n')
	f.write('				</kineticLaw> \n')
	f.write('			</reaction> \n')
	

def write_tail(f):
	f.write('		</listOfReactions> \n')
	f.write('	</model> \n')
	f.write('</sbml>  \n')


def initialize_DNA_amount(num, pos):
    amount=[1]*(num+1);
    for i in pos:
        amount[i]=0;
    return amount;

def initialize_NCoil_amount(num, pos):
    amount=[0]*(num+1);
    for i in pos:
        amount[i]=NCoil_initial;
    return amount;

def initialize_unbind_amount(num, pos):
	amount=[0]*(num+1);
	for i in pos:
		amount[i]=1;
	return amount;

def write_update_step(f, index, k, name_dict, elongation_rate):
	k = k-1; 
	Reactants=[name_dict['RNAP'][k]];
	Products=[name_dict['RNAP'][k], name_dict['step'][k]];
	Propensity_sp=Reactants;
	Propensity_para=elongation_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_pre_elongation(f, index, k, name_dict, Inf):
	k = k-1;
	Reactants=[name_dict['RNAP'][k], name_dict['step'][k]];
	Products=[name_dict['RNAP_tmp'][k]];
	Propensity_sp=Reactants;
	Propensity_para=Inf;
	Reactants_stoi=[1, 60];
	Products_stoi=[1];
	if (k==Z3_index-1):
		Products=[name_dict['RNAP_tmp'][k], 'S2'];
		Products_stoi=[1, 1];
	write_reaction3(f, index, Reactants,Reactants_stoi, Products, Products_stoi, Propensity_sp, Propensity_para);

def write_elongation1(f, index, k, name_dict, Inf): #translocation
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['DNA'][k+1], name_dict['RNAP_tmp'][k], name_dict['DNA'][k+2]];
	Products=[name_dict['DNA'][k], name_dict['RNAP'][k+1], name_dict['DNA'][k+2], name_dict['mRNA'][k]];
	Propensity_para=Inf;
	write_reaction_RNAP_translocation(f, index, Reactants, Products, Propensity_para);


def write_elongation1_2(f, index, k, name_dict, Inf):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['DNA'][k-1], name_dict['RNAP_tmp'][k], name_dict['DNA'][k-2]];
	Products=[name_dict['DNA'][k], name_dict['RNAP'][k-1], name_dict['DNA'][k-2]];
	Propensity_para=Inf;
	write_reaction_RNAP_translocation(f, index, Reactants, Products, Propensity_para);


def write_elongation3(f, index, k, name_dict, Inf): #push all forward
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['Twist'][k]];
	Products=[name_dict['Twist'][k+1]];
	Propensity_sp=[name_dict['RNAP'][k], name_dict['Twist'][k]];
	Propensity_para=Inf*Inf;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_elongation4(f, index, k, name_dict, Inf): 
	k = k-1; 
	Reactants=[name_dict['Twist'][k]];
	Products=[name_dict['Twist'][k+1]];
	Propensity_sp=[name_dict['RNAP_stall'][k], name_dict['Twist'][k]];
	Propensity_para=Inf*Inf;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

'''
def write_elongation4(f, index, k, name_dict, Inf): #no backward
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['Twist'][k]];
	Products=[name_dict['Twist'][k-1]];
	Propensity_sp=[name_dict['RNAP'][k], name_dict['Twist'][k]];
	Propensity_para=0;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);
'''

def write_elongation3_2(f, index, k, name_dict, Inf):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['Twist'][k]];
	Products=[name_dict['Twist'][k-1]];
	Propensity_sp=[name_dict['RNAP'][k], name_dict['Twist'][k]];
	Propensity_para=Inf;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_elongation4_2(f, index, k, name_dict, Inf):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['Twist'][k]];
	Products=[name_dict['Twist'][k+1]];
	Propensity_sp=[name_dict['RNAP'][k], name_dict['Twist'][k]];
	Propensity_para=0;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_elongation_stall(f, index, k, name_dict, lk0, Inf, PCoil_threshold, NCoil_threshold):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['RNAP'][k], name_dict['Twist'][k-1], name_dict['Twist'][k+1]];
	Products=[name_dict['RNAP_stall'][k], name_dict['Twist'][k-1], name_dict['Twist'][k+1]];
	write_reaction_stall(f, index, Reactants, Products, lk0, Inf, PCoil_threshold, NCoil_threshold);

def write_elongation_resumption(f, index, k, name_dict, lk0, resumeRate, PCoil_threshold, NCoil_threshold):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['RNAP_stall'][k], name_dict['Twist'][k-1], name_dict['Twist'][k+1]];
	Products=[name_dict['RNAP'][k], name_dict['Twist'][k-1], name_dict['Twist'][k+1]];
	write_reaction_resume(f, index, Reactants, Products, lk0, resumeRate, PCoil_threshold, NCoil_threshold);

'''
def write_elongation_stall_2(f, index, k, name_dict, lk0, Inf):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['RNAP'][k], name_dict['Twist'][k+1], name_dict['Twist'][k-1]];
	Products=[name_dict['RNAP_stall'][k], name_dict['Twist'][k+1], name_dict['Twist'][k-1]];
	write_reaction_stall(f, index, Reactants, Products, lk0, Inf);

def write_elongation_resumption_2(f, index, k, name_dict, lk0, resumeRate):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['RNAP_stall'][k], name_dict['Twist'][k+1], name_dict['Twist'][k-1]];
	Products=[name_dict['RNAP'][k], name_dict['Twist'][k+1], name_dict['Twist'][k-1]];
	write_reaction_resume(f, index, Reactants, Products, lk0, resumeRate);
'''

def write_premature_termination(f, index, k, name_dict, premature_rate):
	k = k-1; #python uses 0-index system
	Reactants=[name_dict['RNAP_stall'][k]];
	Products=[name_dict['DNA'][k], 'S_pre'];
	Propensity_sp=Reactants;
	Propensity_para=premature_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_elongation_stall_diffusion_D2U(f, index, k, name_dict, stall_diffusion_rate):
	k = k-1; 
	Reactants=[name_dict['RNAP_stall'][k], name_dict['Twist'][k+1]];
	Products=[name_dict['RNAP_stall'][k], name_dict['Twist'][k-1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);
#	write_reaction_stall_diffusion(f, index, Reactants, Products, stall_diffusion_rate1);

def write_elongation_stall_diffusion_U2D(f, index, k, name_dict, stall_diffusion_rate):
	k = k-1;
	Reactants=[name_dict['RNAP_stall'][k], name_dict['Twist'][k-1]];
	Products=[name_dict['RNAP_stall'][k], name_dict['Twist'][k+1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);
#	write_reaction_stall_diffusion(f, index, Reactants, Products, stall_diffusion_rate2);

def write_elongation_translocation_diffusion_D2U(f, index, k, name_dict, stall_diffusion_rate):
	k = k-1;
	Reactants=[name_dict['RNAP'][k], name_dict['Twist'][k+1]];
	Products=[name_dict['RNAP'][k], name_dict['Twist'][k-1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_elongation_translocation_diffusion_U2D(f, index, k, name_dict, stall_diffusion_rate):
	k = k-1;
	Reactants=[name_dict['RNAP'][k], name_dict['Twist'][k-1]];
	Products=[name_dict['RNAP'][k], name_dict['Twist'][k+1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_elongation_translocation_diffusion_D2U_tmp(f, index, k, name_dict, stall_diffusion_rate):
	k = k-1;
	Reactants=[name_dict['RNAP_tmp'][k], name_dict['Twist'][k+1]];
	Products=[name_dict['RNAP_tmp'][k], name_dict['Twist'][k-1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_elongation_translocation_diffusion_U2D_tmp(f, index, k, name_dict, stall_diffusion_rate):
	k = k-1;
	Reactants=[name_dict['RNAP_tmp'][k], name_dict['Twist'][k-1]];
	Products=[name_dict['RNAP_tmp'][k], name_dict['Twist'][k+1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_arrested_promoter_diffusion_D2U(f, index, promoter_index, name_dict, stall_diffusion_rate):
	k = promoter_index-1;
	Reactants=['S4', name_dict['Twist'][k+1]];
	Products=['S4', name_dict['Twist'][k-1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_arrested_promoter_diffusion_U2D(f, index, promoter_index, name_dict, stall_diffusion_rate):
	k = promoter_index-1;
	Reactants=['S4', name_dict['Twist'][k-1]];
	Products=['S4', name_dict['Twist'][k+1]];
	Propensity_sp=Reactants;
	Propensity_para=stall_diffusion_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_initiation(f, index, name_dict, promoter_index, lk0, Boundary0, Boundary1, initiation_rate, mRNA_index):
	Reactants=[name_dict['DNA'][promoter_index-1], name_dict['Twist'][promoter_index-1], name_dict['DNA'][promoter_index]];
	Products=[name_dict['RNAP'][promoter_index-1], name_dict['Twist'][promoter_index-1], name_dict['DNA'][promoter_index], name_dict['S1'][mRNA_index]];
	initiationRate=initiation_rate;
	basalLevel = 0.001/initiation_rate;
	write_reaction_initiation(f, index, Reactants, Products, Boundary0, Boundary1, initiationRate, lk0, basalLevel);

def write_promoter_inhibition(f, index, name_dict, signal_index, promoter_index, Inf):
	Reactants=[name_dict['RNAP'][signal_index-1], name_dict['DNA'][promoter_index-1], name_dict['DNA'][promoter_index]];
	Products=[name_dict['RNAP'][signal_index-1], 'S4', name_dict['DNA'][promoter_index]];
	Propensity_para=Inf;
	write_PI_reaction(f, index, Reactants, Products, Inf);

def write_initiation_2(f, index, name_dict, promoter_index, lk0, initiation_rate, mRNA_index):
	Reactants=[name_dict['DNA'][promoter_index-1], name_dict['Twist'][promoter_index-1], name_dict['DNA'][promoter_index-2]];
	Products=[name_dict['RNAP'][promoter_index-1], name_dict['Twist'][promoter_index-1], name_dict['DNA'][promoter_index-2], name_dict['S1'][mRNA_index]];
	initiationRate=initiation_rate;
	Boundary0 = 0;
	Boundary1 = -0.10;
	basalLevel = 0.005;
	write_reaction_initiation(f, index, Reactants, Products, Boundary0, Boundary1, initiationRate, lk0, basalLevel);

def write_termination(f, index, name_dict, terminator_index, mRNA_index, Inf):
	Reactants=[name_dict['RNAP'][terminator_index-1]];
	Products=[name_dict['DNA'][terminator_index-1]];
	Propensity_sp=Reactants;
	Propensity_para=Inf;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

'''
def write_termination_2(f, index, name_dict, terminator_index, mRNA_index, Inf):
	Reactants=[name_dict['RNAP'][terminator_index-1]];
	Products=[name_dict['DNA'][terminator_index-1], name_dict['mRNA'][mRNA_index], name_dict['S2'][mRNA_index]];
	Propensity_sp=Reactants;
	Propensity_para=Inf;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);
'''
#############
'''
def write_mRNA_degradation(f, index, name_dict, mRNA_index, mRNA_degradation_rate):
	Reactants=[name_dict['mRNA'][mRNA_index]];
	Products="NA";
	Propensity_sp=Reactants;
	Propensity_para=mRNA_degradation_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);
'''

def write_mRNA_degradation_initiation(f, index, name_dict, promoter_index, mRNA_degradation_rate):
	Reactants=[name_dict['mRNA'][promoter_index-1]];
	Products=[name_dict['mRNA_degr'][promoter_index]];
	Propensity_sp=Reactants;
	Propensity_para=mRNA_degradation_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_mRNA_degradation_initiation_pre(f, index, name_dict, promoter_index, mRNA_degradation_rate):
	Reactants=[name_dict['mRNA_pre'][promoter_index-1]];
	Products="NA";
	Propensity_sp=Reactants;
	Propensity_para=mRNA_degradation_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_mRNA_degradation_elongation(f, index, name_dict, k, RNase_elongation_rate):
	k = k-1;
	Reactants=[name_dict['mRNA'][k], name_dict['mRNA_degr'][k], name_dict['mRNA_pre'][k]];
	Products=[name_dict['mRNA_degr'][k+1], name_dict['mRNA_pre'][k]];
	write_reaction_mRNA_degradation_elongation(f, index, Reactants, Products, RNase_elongation_rate);

def write_mRNA_degradation_elongation_pre(f, index, name_dict, k, RNase_elongation_pre_rate):
	k = k-1;
	Reactants=[name_dict['mRNA_pre'][k], name_dict['mRNA_degr'][k]];
	Products="NA";
	Propensity_sp=Reactants;
	Propensity_para=RNase_elongation_pre_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_mRNA_degradation_termination(f, index, name_dict, terminator_index, RNase_elongation_rate):
	Reactants=[name_dict['mRNA'][terminator_index-1], name_dict['mRNA_degr'][terminator_index-1]];
	Products="NA";
	Propensity_sp=Reactants;
	Propensity_para=RNase_elongation_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_translation(f, index, name_dict, terminator_index, mRNA_index, translation_rate):
	Reactants=[name_dict['mRNA'][terminator_index-1]];
	Products=[name_dict['mRNA'][terminator_index-1], name_dict['protein'][mRNA_index]];
	Propensity_sp=Reactants;
	Propensity_para=translation_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_protein_degradation(f, index, name_dict, mRNA_index, protein_degradation_rate):
	Reactants=[name_dict['protein'][mRNA_index]];
	Products="NA";
	Propensity_sp=Reactants;
	Propensity_para=protein_degradation_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

############
def write_diffusion1(f, index, k, name_dict, diffusion_rate):
	k = k-1;
	Reactants=[name_dict['Twist'][k], name_dict['Twist'][k+1], name_dict['DNA'][k+1]];
	Products=[name_dict['Twist'][k+1], name_dict['DNA'][k+1]];
	Products_stoi = [2, 1];
	write_reaction_SC_diffusion(f, index, Reactants, Products, Products_stoi, diffusion_rate);

def write_diffusion2(f, index, k, name_dict, diffusion_rate):
	k = k-1;
	Reactants=[name_dict['Twist'][k], name_dict['Twist'][k-1], name_dict['DNA'][k-1]];
	Products=[name_dict['Twist'][k-1], name_dict['DNA'][k-1]];
	Products_stoi = [2, 1];
	write_reaction_SC_diffusion(f, index, Reactants, Products, Products_stoi, diffusion_rate);


def write_diffusion3(f, index, k, name_dict, diffusion_rate):
	k = k-1;
	Reactants=[name_dict['Twist'][k], name_dict['DNA'][k+1]];
	Products=[name_dict['Twist'][k+1], name_dict['DNA'][k+1]];
	Propensity_para=diffusion_rate;
	Propensity_sp=Reactants;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_diffusion4(f, index, k, name_dict, diffusion_rate):
	k = k-1;
	Reactants=[name_dict['Twist'][k], name_dict['DNA'][k-1]];
	Products=[name_dict['Twist'][k-1], name_dict['DNA'][k-1]];
	Propensity_para=diffusion_rate;
	Propensity_sp=Reactants;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);


def write_diffusion5(f, index, k, name_dict, diffusion_rate):
	Reactants=[name_dict['Twist'][0], name_dict['Twist'][num-1], name_dict['DNA'][num-1]];
	Products=[name_dict['Twist'][num-1], name_dict['DNA'][num-1]];
	Products_stoi = [2, 1];
	write_reaction_SC_diffusion(f, index, Reactants, Products, Products_stoi, diffusion_rate);

def write_diffusion6(f, index, k, name_dict, rate2):
	Reactants=[name_dict['Twist'][num-1], name_dict['Twist'][0], name_dict['DNA'][0]];
	Products=[name_dict['Twist'][0], name_dict['DNA'][0]];
	Products_stoi = [2, 1];
	write_reaction_SC_diffusion(f, index, Reactants, Products, Products_stoi, diffusion_rate);

# birth
def write_diffusion7(f, index, k, name_dict, rate1):
	k = k-1;
	Reactants=[name_dict['DNA'][k]];
	Products=[name_dict['DNA'][k], name_dict['Twist'][k]];
	Propensity_para=rate1;
	Propensity_sp=Reactants;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

#death 
def write_diffusion8(f, index, k, name_dict, rate2):
	k = k-1;
	Reactants=[name_dict['Twist'][k]];
	Products="NA";
	Propensity_para=rate2;
	Propensity_sp=Reactants;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

### gyrase bind and unbind
def write_gyrase_binding1(f, index, k, name_dict, gyrase_binding_rate):
	k = k-1;
	Reactants=[name_dict['Gyrase_unbind'][k], name_dict['DNA'][k]];
	Products=[name_dict['Gyrase_bind'][k], name_dict['DNA'][k]];
	Propensity_sp=Reactants;
	Propensity_para=gyrase_binding_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_gyrase_binding2(f, index, k, name_dict, NCoil_producing_rate, lk0):
	k = k-1;
	Reactants=[name_dict['Gyrase_bind'][k], name_dict['Twist'][k]];
	Products=[name_dict['Gyrase_bind'][k]];
	sigma_t = 0.0145;
#	sigma_t = 0.01;
	epsilon = 0.0122;
	catalyticRate = NCoil_producing_rate;
	write_reaction_gyrase_catalysis(f, index, Reactants, Products, sigma_t, epsilon, catalyticRate, lk0);

def write_gyrase_unbind(f, index, k, name_dict, gyrase_unbinding_rate):
	k = k-1;
	Reactants=[name_dict['Gyrase_bind'][k]];
	Products=[name_dict['Gyrase_unbind'][k]];
	Propensity_sp=Reactants;
	Propensity_para=gyrase_unbinding_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

### topoI binding and unbind
def write_topoI_binding1(f, index, k, name_dict, topoI_binding_rate):
	k = k-1;
	Reactants=[name_dict['TopoI_unbind'][k], name_dict['DNA'][k]];
	Products=[name_dict['TopoI_bind'][k], name_dict['DNA'][k]];
	Propensity_sp=Reactants;
	Propensity_para=topoI_binding_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_topoI_binding2(f, index, k, name_dict, PCoil_producing_rate, lk0):
	k = k-1;
	Reactants=[name_dict['TopoI_bind'][k], name_dict['Twist'][k]];
	Products=[name_dict['TopoI_bind'][k], name_dict['Twist'][k]];
	Products_stoi = [1, 2];
	catalyticRate = PCoil_producing_rate;
	write_reaction_TopoI_catalysis(f, index, Reactants, Products, Products_stoi, catalyticRate, lk0);

def write_topoI_unbind(f, index, k, name_dict, topoI_unbinding_rate):
	k = k-1;
	Reactants=[name_dict['TopoI_bind'][k]];
	Products=[name_dict['TopoI_unbind'][k]];
	Propensity_sp=Reactants;
	Propensity_para=topoI_unbinding_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

### loop and unloop
def write_loop1(f, index, loop_rate): # from unloop to loop
	Reactants=['unloop_state'];
	Products=['loop1', 'loop2', 'loop_state'];
	Propensity_sp=Reactants;
	Propensity_para=loop_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_loop2(f, index, loop_sp, DNA_sp, Inf):
	Reactants=[loop_sp, DNA_sp]; #'loop1', 'DNA11'; 'loop2', 'DNA38'
	Products="NA"; # 'DNA11' and 'DNA38' is occupied
	Propensity_sp=Reactants;
	Propensity_para=Inf;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_unloop1(f, index, unloop_rate): # from loop to unloop
	Reactants=['loop_state'];
	Products=['unloop1', 'unloop2', 'unloop_state'];
	Propensity_sp=Reactants;
	Propensity_para=unloop_rate;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);

def write_unloop2(f, index, unloop_sp, DNA_sp, Inf):
	Reactants=[unloop_sp];
	Products=[DNA_sp]; # 'DNA11' and 'DNA38' is free
	Propensity_sp=Reactants;
	Propensity_para=Inf;
	write_reaction(f, index, Reactants, Products, Propensity_sp, Propensity_para);


###########################################
# rate constant
#basic 
lk0 = 60;

#initiation_rate=float(sys.argv[1]);
Inf=10000000000000000000000000;
elongation_rate=60;
diffusion_rate=500;
#diffusion_rate=500;
mRNA_degradation_rate=0.0067;
RNase_elongation_rate=Inf;
RNase_elongation_pre_rate=Inf*Inf;
translation_rate=0.03;
protein_degradation_rate=0.002;
resumeRate=Inf;

# loop
loop_rate=0;
unloop_rate=0;
#loop_rate=0.00833;  # every 2 min
#unloop_rate=0.00167;# every 10 min

# Gyrase
gyrase_binding_rate=0.0018;    # every 1s
#gyrase_binding_rate=0;
gyrase_unbinding_rate=0.5;
NCoil_producing_rate=20;

# topoI
topoI_binding_rate=0.0018;
#topoI_binding_rate=0;
topoI_unbinding_rate=0.5;
PCoil_producing_rate=10;

# NCoil
relaxation_rate1=50;
stall_diffusion_rate=0.2;
#relaxation_rate1=float(sys.argv[1]); #1,5,10,50
#stall_diffusion_rate=float(sys.argv[2]);

NCoil_initial=56;
rate1=56*relaxation_rate1;
rate2=relaxation_rate1;

rate11=rate1;
rate22=rate2;

#stall
initiation_rate=float(sys.argv[1]);
initiation_list = [initiation_rate];
SC_threshold=0.6;
premature_rate=0;
#SC_threshold=float(sys.argv[1]);
#premature_rate=float(sys.argv[2]);
NCoil_threshold=-SC_threshold;
PCoil_threshold=SC_threshold;

###########################################
# global config:
'''
num=88;
barrier=[];
promoter_list=[2];
terminator_list=[87];
signal_index=47;
Z3_index=52;
'''

num=149;
barrier=[];
promoter_list=[30];
terminator_list=[115];
signal_index=75;
Z3_index=80;


'''
free_NCoil=range(1,barrier[0]);
free_NCoil.extend(range(barrier[0]+1,barrier[1]));
free_NCoil.extend(range(barrier[1]+1,num+1));
gyrase_binding_site=range(1, num+1);
topoI_binding_site=range(1, num+1);
'''

free_NCoil=range(1, num+1);
gyrase_binding_site=range(1, num+1);
topoI_binding_site=range(1, num+1);

name_dict={};
for i in ['DNA', 'RNAP','PCoil','NCoil', 'RNAP_stall', 'Gyrase_unbind', 'Gyrase_bind', 'TopoI_unbind', 'TopoI_bind', 'RNAP_tmp', 'step', 'Twist', 'mRNA', 'mRNA_degr', 'mRNA_pre']:
	name_dict[i] = [i+str(j) for j in range(1,num+1)];

mRNA_index=[0];
name_dict['protein'] = ['protein1'];
name_dict['S1'] = ['S1'];
name_dict['S2'] = ['S2'];

Boundary0_list = [0];
Boundary1_list = [-0.06];
###########################################
# write files
filename='expl_sens_500.'+str(sys.argv[1])+'.sbml';
f = open(filename, 'w')
wirte_head(f);
DNA_amount=initialize_DNA_amount(num, barrier);
for i in range(0,num):					#0-1
	name=name_dict['DNA'][i];
	write_sp(f, name, DNA_amount[i+1]);

for j in ['RNAP', 'RNAP_stall']:		#1-2, 2-3
	for i in range(0,num):
		name=name_dict[j][i];
		write_sp(f, name, 0);


NCoil_amount=initialize_NCoil_amount(num, free_NCoil);
#NCoil_amount[1]=0;
#NCoil_amount[num]=0;

for i in range(0,num):					#3-4
	name=name_dict['Twist'][i];
	write_sp(f, name, NCoil_amount[i+1]);

# gyrase related: 
unbind_amount=initialize_unbind_amount(num, gyrase_binding_site);
for i in range(0,num):					#4-5
	name=name_dict['Gyrase_unbind'][i];
	write_sp(f, name, unbind_amount[i+1]);
for i in range(0,num):					#5-6
	name=name_dict['Gyrase_bind'][i];
	write_sp(f, name, 0);
# topoI related:
unbind_amount=initialize_unbind_amount(num, topoI_binding_site);
for i in range(0,num):					#6-7
	name=name_dict['TopoI_unbind'][i];
	write_sp(f, name, unbind_amount[i+1]);
for i in range(0,num):					#7-8
	name=name_dict['TopoI_bind'][i];
	write_sp(f, name, 0);

#RNAP_tmp
for i in range(0,num):					#8-9
	name=name_dict['RNAP_tmp'][i];
	write_sp(f, name, 0);

for i in range(0,num):					#9-10
	name=name_dict['step'][i];
	write_sp(f, name, 0);

# mRNA & mRNA_degradation
for i in range(0,num):					#10-11
	name=name_dict['mRNA'][i];
	write_sp(f, name, 0);

for i in range(0,num):					#11-12
	name=name_dict['mRNA_degr'][i];
	write_sp(f, name, 0);
for i in range(0,num): 
	name=name_dict['mRNA_pre'][i];
	write_sp(f, name, 0);

# protein & probe stuff
for i in range(0,len(name_dict['protein'])):	#12+0
	name=name_dict['protein'][i];
	write_sp(f, name, 0);

for i in range(0,len(name_dict['protein'])):	#12+1
	name=name_dict['S1'][i];
	write_sp(f, name, 0);

for i in range(0,len(name_dict['protein'])):	#12+2
	name=name_dict['S2'][i];
	write_sp(f, name, 0);

write_sp(f, 'S4_pre', 1);
write_sp(f, 'S4', 0);
write_sp(f, 'S_pre', 0);
# wirte others
#write_sp(f, 'Z1', 0);
#write_sp(f, 'Z2', 0);
#write_sp(f, 'Gyrase', gyrase_num);
#write_sp(f, 'TopoI', topoI_num);
#write_sp(f, 'RNAP', RNAP_num);
# write loop and unloop
write_sp(f, 'loop_state', 1);
write_sp(f, 'unloop_state', 0);
write_sp(f, 'loop1', 0);
write_sp(f, 'loop2', 0);
write_sp(f, 'unloop1', 0);
write_sp(f, 'unloop2', 0);
f.write('		</listOfSpecies> \n')
f.write('		<listOfReactions> \n')

index=1;
# transcription
#'''
for i in range(0,len(promoter_list)):
	promoter_i = promoter_list[i];
	terminator_i = terminator_list[i];
	initiation_rate_i = initiation_list[i];
	Boundary0_i = Boundary0_list[i];
	Boundary1_i = Boundary1_list[i];
	mRNA_i = mRNA_index[i];
	write_initiation(f, index, name_dict, promoter_i, lk0, Boundary0_i, Boundary1_i, initiation_rate_i, mRNA_i); index = index+1;
#	write_promoter_inhibition(f, index, name_dict, signal_index, promoter_i, Inf); index = index+1;
#	write_arrested_promoter_diffusion_D2U(f, index, promoter_i, name_dict, stall_diffusion_rate); index = index+1;
#	write_arrested_promoter_diffusion_U2D(f, index, promoter_i, name_dict, stall_diffusion_rate); index = index+1;
	pos_index = 0;
	for k in range(promoter_i, terminator_i):
		write_update_step(f, index, k, name_dict, elongation_rate); index = index+1;
		write_pre_elongation(f, index, k, name_dict, Inf); index = index+1;
		write_elongation1(f, index, k, name_dict, Inf); index = index+1;
		write_elongation_stall(f, index, k, name_dict, lk0, Inf, PCoil_threshold, NCoil_threshold); index = index+1;
		write_elongation_resumption(f, index, k, name_dict, lk0, resumeRate, PCoil_threshold, NCoil_threshold); index = index+1;
#		write_premature_termination(f, index, k, name_dict, premature_rate); index = index+1;
		write_elongation_stall_diffusion_D2U(f, index, k, name_dict, stall_diffusion_rate); index = index+1;
		write_elongation_stall_diffusion_U2D(f, index, k, name_dict, stall_diffusion_rate); index = index+1;
		write_elongation_translocation_diffusion_D2U(f, index, k, name_dict, stall_diffusion_rate); index = index+1;
		write_elongation_translocation_diffusion_U2D(f, index, k, name_dict, stall_diffusion_rate); index = index+1;
		write_elongation_translocation_diffusion_D2U_tmp(f, index, k, name_dict, stall_diffusion_rate); index = index+1;
		write_elongation_translocation_diffusion_U2D_tmp(f, index, k, name_dict, stall_diffusion_rate); index = index+1;
		pos_index += 1;
	for k in range(promoter_i, terminator_i+1):
		write_elongation3(f, index, k, name_dict, Inf); index = index+1;
		write_elongation4(f, index, k, name_dict, Inf); index = index+1;
	write_termination(f, index, name_dict, terminator_i, mRNA_i, Inf); index = index+1;
	write_mRNA_degradation_initiation(f, index, name_dict, promoter_i, mRNA_degradation_rate); index = index+1;
	write_mRNA_degradation_initiation_pre(f, index, name_dict, promoter_i, mRNA_degradation_rate); index = index+1;
	for k in range(promoter_i+1, terminator_i-1):
		write_mRNA_degradation_elongation(f, index, name_dict, k, RNase_elongation_rate); index = index+1;
		write_mRNA_degradation_elongation_pre(f, index, name_dict, k, RNase_elongation_pre_rate); index = index+1;
	write_mRNA_degradation_termination(f, index, name_dict, terminator_i-1, RNase_elongation_rate); index = index+1;
	write_mRNA_degradation_elongation_pre(f, index, name_dict, terminator_i-1, RNase_elongation_pre_rate); index = index+1;
	write_translation(f, index, name_dict, terminator_i-1, mRNA_i, translation_rate); index = index+1;
	write_protein_degradation(f, index, name_dict, mRNA_i, protein_degradation_rate);
#'''
'''
for i in range(0,len(promoter_index2)):
	promoter_i = promoter_index2[i];
	terminator_i = terminator_index2[i];
	mRNA_i = mRNA_index2[i];
	write_initiation_2(f, index, name_dict, promoter_i, lk0, initiation_rate, mRNA_i); index = index+1;
	for k in reversed(range(terminator_i+1, promoter_i+1)):
		write_update_step(f, index, k, name_dict, elongation_rate); index = index+1;
		write_pre_elongation(f, index, k, name_dict, Inf); index = index+1;
		write_elongation1_2(f, index, k, name_dict, Inf); index = index+1;
		write_elongation_stall_2(f, index, k, name_dict, lk0, Inf); index = index+1;
		write_elongation_stall_diffusion_D2U(f, index, k, name_dict, stall_diffusion_rate1); index = index+1;
		write_elongation_stall_diffusion_U2D(f, index, k, name_dict, stall_diffusion_rate2); index = index+1;
		write_elongation_resumption_2(f, index, k, name_dict, lk0, resumeRate); index = index+1;
	for k in reversed(range(terminator_i, promoter_i+1)):
		write_elongation3_2(f, index, k, name_dict, Inf); index = index+1;
		write_elongation4_2(f, index, k, name_dict, Inf); index = index+1;
	write_termination_2(f, index, name_dict, terminator_i, mRNA_i, Inf); index = index+1;
	write_mRNA_degradation(f, index, name_dict, mRNA_i, mRNA_degradation_rate); index = index+1;
	write_translation(f, index, name_dict, mRNA_i, translation_rate); index = index+1;
	write_protein_degradation(f, index, name_dict, mRNA_i, protein_degradation_rate);
'''
# diffusion
for k in range(1,num):
	write_diffusion3(f, index, k, name_dict, diffusion_rate); index = index+1;
for k in range(2,num+1):
	write_diffusion4(f, index, k, name_dict, diffusion_rate); index = index+1;
#write_diffusion5(f, index, k, name_dict, diffusion_rate); index = index+1;
#write_diffusion6(f, index, k, name_dict, diffusion_rate); index = index+1;

write_diffusion7(f, index, 1, name_dict, rate1); index = index+1;
write_diffusion8(f, index, 1, name_dict, rate2); index = index+1;

write_diffusion7(f, index, num, name_dict, rate11); index = index+1;
write_diffusion8(f, index, num, name_dict, rate22); index = index+1;

# topomerase
for k in gyrase_binding_site:
	write_gyrase_binding1(f, index, k, name_dict, gyrase_binding_rate); index = index+1;
	write_gyrase_binding2(f, index, k, name_dict, NCoil_producing_rate, lk0); index = index+1;
	write_gyrase_unbind(f, index, k, name_dict, gyrase_unbinding_rate); index = index+1;
for k in topoI_binding_site:
	write_topoI_binding1(f, index, k, name_dict, topoI_binding_rate); index = index+1;
	write_topoI_binding2(f, index, k, name_dict, PCoil_producing_rate, lk0); index = index+1;
	write_topoI_unbind(f, index, k, name_dict, topoI_unbinding_rate); index = index+1;
# loop and unloop

'''
B1 = 'DNA'+str(barrier[0])
B2 = 'DNA'+str(barrier[1])
write_loop1(f, index, loop_rate); index = index+1;
write_loop2(f, index, 'loop1', B1, Inf); index = index+1;
write_loop2(f, index, 'loop2', B2, Inf); index = index+1;
write_unloop1(f, index, unloop_rate); index = index+1;
write_unloop2(f, index, 'unloop1', B1, Inf); index = index+1;
write_unloop2(f, index, 'unloop2', B2, Inf); index = index+1;
'''
write_tail(f);

