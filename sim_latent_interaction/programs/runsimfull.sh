#!/bin/bash
#$ -cwd
#$ -o  /u/project/cenders/remusmit/sim_latent_interaction/joboutputs/$JOB_ID.txt
#$ -j y
. /u/local/etc/profile.d/sge.sh #qsub command wasn't found until i added this line
RUNONCLUSTER=0

for CATLOOP in {0..1}; do
	for PROBLOOP in {0..2}; do #3
		for RSQLOOP in {0..2}; do #3
			for NLOOP in {0..8}; do #9
				for LOADLOOP in {0..1}; do #2
					for ITEMLOOP in {0..1}; do #2



						if [ ${RUNONCLUSTER} = 0 ]
						then
							PROGDIR=$(dirname "$BASH_SOURCE")
							DIRNAME="$(dirname "$PROGDIR")"
							export RUNONCLUSTER CATLOOP PROBLOOP RSQLOOP NLOOP LOADLOOP ITEMLOOP
							/bin/bash ${PROGDIR}/runsim.sh
						else
							
							qsub -v RUNONCLUSTER=${RUNONCLUSTER},CATLOOP=${CATLOOP},PROBLOOP=${PROBLOOP},RSQLOOP=${RSQLOOP},NLOOP=${NLOOP},LOADLOOP=${LOADLOOP},ITEMLOOP=${ITEMLOOP} < runsim.sh 
						fi
					done
	    		done
			done
		done
	done
done

