#!/bin/bash
#$ -cwd
#$ -o  /u/project/cenders/remusmit/sim_latent_interaction/joboutputs/$JOB_ID.txt
#$ -j y
#  Resources requested
#$ -l h_data=10G,h_rt=3:00:00
#$ -t 1-10:1


NUMREPS=10

cat=(2 3)
groupprob=(1 2 3)
rsq=(0 .03 .07)
sample=(100  150  200  250  300  350  400  500 1000)
loading=(.5 .8)
nitem=(6 12)

# cluster specifications
INTERACTIVE=1
MULTIPLE=1

if [ ${INTERACTIVE} = 1 ]
	then
		RUNONCLUSTER=1
		CATLOOP=0
		PROBLOOP=0
		RSQLOOP=0
		NLOOP=0
		LOADLOOP=0
	    ITEMLOOP=1

	fi

if [ ${RUNONCLUSTER} = 0 ]
then
	# paths
	DIRNAME=~/Documents/GitHub/latent_interaction/sim_latent_interaction
	PROGDIR=${DIRNAME}/programs
	MISCDIR=${DIRNAME}/misc
	OUTFDIR=${DIRNAME}/joboutputs
	RESDIR=${DIRNAME}/results

	# application paths
	MPATH=/applications/mplus/mplus
	RPATH=R
	BPATH=/Applications/Blimp/blimp
	SGE_TASK_ID=1

else	 
	# PATHS
	DIRNAME=/u/project/cenders/remusmit/sim_latent_interaction
	PROGDIR=${DIRNAME}/programs
	MISCDIR=${DIRNAME}/misc
	RESDIR=${DIRNAME}/results
	OUTFDIR=${DIRNAME}/joboutputs
	BPATH=blimp
	RPATH=R

	# LOAD MODULES
	. /u/local/Modules/default/init/modules.sh
	source /u/local/Modules/default/init/modules.sh #IDRE support said to add this line
	module use /u/project/cenders/apps/modulefiles
	module load gcc/10.2.0
	module load blimp/3.2.3
	module load nlopt
	module load R/4.4.0

fi
# housekeeping: clean out folders before starting sim
# rm ${OUTFDIR}/*.*
# rm ${RESDIR}/*.*


# REP RANGE
REPFIRST=${SGE_TASK_ID}
REPLAST=$((${REPFIRST} + ${MULTIPLE}))

# print start time to determine runtime
#echo "conditions: gprobs = ${PROBLOOP}, RSQ = ${RSQLOOP}, sample = ${NLOOP}, loading = ${LOADLOOP}" #>> ${MISCDIR}/time_log.txt
#echo "start time: " ` date ` #>> ${MISCDIR}/time_log.txt



for (( REPLOOP = ${REPFIRST}; REPLOOP < ${REPLAST}; REPLOOP++ )); do

	FILENAME=cat${CATLOOP}prob${PROBLOOP}rsq${RSQLOOP}N${NLOOP}load${LOADLOOP}item${ITEMLOOP}rep${REPLOOP}

	SEEDLINECAT=$((${CATLOOP} * ${#groupprob[@]} * ${#rsq[@]} * ${#sample[@]} * ${#loading[@]} * ${#nitem[@]} + ${REPLOOP} + 1))
	SEEDLINEPR=$((${PROBLOOP} * ${#rsq[@]} * ${#sample[@]} * ${#loading[@]} * ${#nitem[@]} + ${REPLOOP}))
	SEEDLINER=$((${RSQLOOP} * ${#sample[@]} * ${#loading[@]} * ${#nitem[@]} + ${REPLOOP}))
	SEEDLINEN=$((${NLOOP} * ${#loading[@]} * ${#nitem[@]} + ${REPLOOP}))
	SEEDLINEL=$((${LOADLOOP} * ${#nitem[@]} + ${REPLOOP}))
	SEEDLINENI=$((${ITEMLOOP} + ${REPLOOP}))
	SEEDLINENUM=$((${SEEDLINECAT} + ${SEEDLINEPR} + ${SEEDLINER} + ${SEEDLINEN} + ${SEEDLINEL} + ${SEEDLINENI}))


	SEED=$(sed -n "${SEEDLINENUM}p" "${MISCDIR}/seedlist.dat")

	${RPATH} --no-save --slave --args ${RUNONCLUSTER} ${DIRNAME} ${FILENAME} ${cat[CATLOOP]} ${groupprob[PROBLOOP]} ${rsq[RSQLOOP]} ${sample[NLOOP]} ${loading[LOADLOOP]} ${nitem[ITEMLOOP]} ${REPLOOP} ${SEED} < ${PROGDIR}/simulation_one_rep.R

	#echo "${RUNONCLUSTER} ${DIRNAME} ${cat[CATLOOP]} ${groupprob[PROBLOOP]} ${rsq[RSQLOOP]} ${sample[NLOOP]} ${loading[LOADLOOP]} ${nitem[ITEMLOOP]} ${REPLOOP} ${SEED}"

done   

# print end time to determine runtime
#echo "end time: " ` date ` #>> ${MISCDIR}/time_log.txt


