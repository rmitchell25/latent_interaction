#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o  /u/project/cenders/remusmit/sim_latent_interaction/joboutputs/$JOB_ID.txt
#$ -j y



NUMREPS=2000

cat=1
groupprob=2
rsq=2
sample=9
loading=1
nitem=1

FOLDERPATH=/u/home/r/remusmit/project-cenders/sim_latent_interaction
SCRATCHPATH=/u/scratch/r/remusmit
SEEDLIST=${FOLDERPATH}/misc/seedlist.dat

for (( CONDCAT=0; CONDCAT<2; CONDCAT++ )); do
	for (( CONDGROUP=0; CONDGROUP<3; CONDGROUP++ )); do
		for (( CONDRSQ=0; CONDRSQ<3; CONDRSQ++ )); do
			for (( CONDSAMPLE=0; CONDSAMPLE<9; CONDSAMPLE++ )); do
				for (( CONDLOAD=0; CONDLOAD<2; CONDLOAD++ )); do
					for (( CONDNITEM=0; CONDNITEM<2; CONDNITEM++ )); do
						for (( REP=1; REP<2001; REP++ )); do

							FILENAME=cat${CONDCAT}prob${CONDGROUP}rsq${CONDRSQ}N${CONDSAMPLE}load${CONDLOAD}item${CONDNITEM}rep${REP}.dat

							echo "${FILENAME}" 
							

							if [ -f ${SCRATCHPATH}/${FILENAME} ]; then
								:
							else
								LINENUM=$(( 
									CONDCAT   *(groupprob*rsq*sample*loading*nitem*NUMREPS) + 
									CONDGROUP*(rsq*sample*loading*nitem*NUMREPS) + 
									CONDRSQ  *(sample*loading*nitem*NUMREPS) + 
									CONDSAMPLE*(loading*nitem*NUMREPS) + 
									CONDLOAD *(nitem*NUMREPS) + 
									CONDNITEM*(NUMREPS) + 
									REP ))

								SEED=$(sed -n "${LINENUM}p" "${SEEDLIST}")

								#echo "${FILENAME}"

								echo "cat${CONDCAT}prob${CONDGROUP}rsq${CONDRSQ}N${CONDSAMPLE}load${CONDLOAD}item${CONDNITEM}rep${REP}" >> ${FOLDERPATH}/misc/missingreplist.txt
								echo "${SEED}" >> ${FOLDERPATH}/misc/missingseedlist.txt
							fi

						done
					done
				done
			done
		done
	done
done
