#!/bin/zsh
# This script submits all jobs to naf, which are defined in dataset

version="04"
datasets=(

# QCD
 /pnfs/desy.de/cms/tier2/store/user/jschulz/nTuples/GJets_40_100_V09/
 /pnfs/desy.de/cms/tier2/store/user/jschulz/nTuples/GJets_100_200_V09/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/GJets_200_400_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/GJets_400_inf_V03/
 /pnfs/desy.de/cms/tier2/store/user/jschulz/nTuples/QCD_100-250_V09/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/QCD_250-500_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/QCD_500-1000_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/QCD_1000-inf_V03/

# EWK
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/TTJets_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/WJets_250_300_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/WJets_300_400_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/WJets_400_inf_V03/

# ISR
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/TTGamma_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/WGamma_50_130_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/WGamma_130_inf_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/ZGammaNuNu_V03/
 /pnfs/desy.de/cms/tier2/store/user/jschulz/nTuples/ZGamma_V02/

## Data
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/PhotonHadA_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/PhotonHadB_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/PhotonHadC_V03/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/PhotonHadD_V03/

 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/PhotonA_V04/ 
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/SinglePhotonB_V04/
 /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/SinglePhotonC_V04/
 /pnfs/desy.de/cms/tier2/store/user/jschulz/nTuples/PhotonParkedD_V10/
 

## Signal/
# /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/T5gg_V04/
# /pnfs/desy.de/cms/tier2/store/user/kiesel/nTuples/T5wg_V04/
)

# settings
outputFolder=/nfs/dust/cms/user/eicker
filesPerJob=20
prefix="dcap://dcache-cms-dcap.desy.de"

for dataset in "${datasets[@]}"; do

    # get folder name as best description for job
    job_name=$(echo $dataset|rev|cut -d'/' -f2|rev)
    # Since root can't handle -, it will be substituted to _
    job_name=$(echo $job_name|sed 's/-/_/g')

    files=( $dataset*root )

    for (( i=0; i<=$(( ${#files[*]} -1 )); i+=$filesPerJob )) {
        thisFiles=${files[@]:$i:$filesPerJob}

        #eval echo ${prefix}\{${thisFiles// /,}\}

        jobPrefix=${job_name}.${version}__${i}

        outputFileName=$outputFolder/${jobPrefix}_tree.root
        
        eval qsub -N ${job_name}_${i} -b y -j y -l os=sld6 -l h_vmem=1000M -l h_rt=11:00:00 -l site=hh \
            `pwd`/submitScriptTemplate.sh $outputFileName ${prefix}\{${thisFiles// /,}\}
        # b: execute a binary or script that will not be transferred to the batch node thus the batch node needs to have access to the given path
        # j: merge stdout and stderr
        # o: log file
    }
    # all jobs
done #dataset
