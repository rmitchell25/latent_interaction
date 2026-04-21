#!/bin/sh
#$ -cwd
#$ -j y
#  Resources requested
#$ -o  /u/project/cenders/remusmit/sim_latent_interaction/joboutputs/$JOB_ID.txt
#$ -l h_data=10G,h_rt=3:00:00
#$ -t 1-75000:1
#$ -l highp

# paths
DIRNAME=/u/home/r/remusmit/project-cenders/sim_latent_interaction
PROGDIR=${DIRNAME}/programs
BPATH=blimp
RPATH=R


# LOAD MODULES
. /u/local/Modules/default/init/modules.sh
source /u/local/Modules/default/init/modules.sh #IDRE support said to add this line
module use /u/project/cenders/apps/modulefiles
module load -f gcc/10.2.0
module load -f blimp
module load nlopt
module load -f R


# define variables
REPLIST=${DIRNAME}/misc/missingreplist.txt
SEEDLIST=${DIRNAME}/misc/missingseedlist.txt
LENGTH=$(wc -l ${REPLIST} | awk '{print $1}')

RUNONCLUSTER=1
cat=(2 3)
groupprob=(1 2 3)
rsq=(0 .03 .07)
sample=(100  150  200  250  300  350  400  500 1000)
loading=(.5 .8)
nitem=(6 12)


# run leftovers
#for (( I =1; I < (${LENGTH}+1); I++ )); do
I=${SGE_TASK_ID}

	FILENAME=$(sed -n "${I}p" "$REPLIST")
	SEED=$(sed -n "${I}p" "$SEEDLIST")
	INFILE=${DIRNAME}/results/${FILENAME}.dat


    cat_level=$(echo "$FILENAME" | sed -E 's/.*cat([0-9]+).*/\1/')
    prob_level=$(echo "$FILENAME" | sed -E 's/.*prob([0-9]+).*/\1/')
    rsq_level=$(echo "$FILENAME" | sed -E 's/.*rsq([0-9]+).*/\1/')
    n_level=$(echo "$FILENAME" | sed -E 's/.*N([0-9]+).*/\1/')
    load_level=$(echo "$FILENAME" | sed -E 's/.*load([0-9]+).*/\1/')
    nitem_level=$(echo "$FILENAME" | sed -E 's/.*item([0-9]+).*/\1/')

	rep_value=$(echo "$FILENAME" | sed -E 's/.*rep([0-9]+).*/\1/') 


    #echo "FILE=${FILENAME}, cat=${cat_level}, groupprob=${prob_level}, effect=${ab_level}, rsq=${rsq_level}, sample=${n_level}, loading=${load_level}, nitem=${nitem_level}, rep=${rep_value}, seed = ${SEED}" >> ${DIRNAME}/misc/check.dat

	${RPATH} --no-save --slave --args ${RUNONCLUSTER} ${DIRNAME} ${FILENAME} ${cat[cat_level]} ${groupprob[prob_level]} ${rsq[rsq_level]} ${sample[n_level]} ${loading[load_level]} ${nitem[nitem_level]}  ${rep_value} ${SEED} < ${PROGDIR}/redo_missing.R
#done

