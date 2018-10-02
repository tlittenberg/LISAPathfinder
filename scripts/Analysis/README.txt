File Structure:

							\Analysis						
								:
			-	-	-	-	-	-	-	-	-	-	-	-	-	-
			:		:						:					:		
			\data	\html_scripts			\plots			\scripts
			:		:						:					:	-	-	-	-	-	-	-	-	-	-	-	
			:		:						-	-	-	-													:
			:	---------								:													:
			:	Directory containing html scripts	---------												:
			:	for viewing individual impacts		A\*SEGMENT_TIME*										:
			:	Dependent on structure in \plots															:
			:																								:
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-						:
	:			:					:					:		:					:						:
\ALL_IMPACTS	\ONLY_IMPACTS	\FLUX_IMPACTS_ALL	\models	\population_models	\population_ratios			:
	:			:					:					:		:					:						:
--------		:				--------				:	--------				:						:
all .pickle						.npy files				:	Pop models, cumul		:						:
files for segments 				for all segments		:	bins					:						:
searched		:				searched			--------					--------					:
				:									micrometeoroid 				Files containing			:
			--------								files reformatted			log likelihood				:
			Contains .pickle files					(from Petr's files)			and ratios from 			:
			for only segments containing			into non-cumul bins			bayesian population 		:
			an impact															inference					:
																											:
																											:
																											:
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
												SCRIPTS
	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
	
	before running these scripts, you will definitely want to check the definition of 
	BASE_DIR : where the \Analysis directory lives



