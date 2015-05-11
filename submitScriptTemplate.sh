#!/bin/bash

echo Running: treeWriter $*

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/desy.de/user/k/kiesel/scratch/singlePhoton/TreeWriterNiklas
#export LD_LIBRARY_PATH=/opt/d-cache/dcap/lib64:/afs/desy.de/products/root/amd64_rhel60/5.34.00/lib:/usr/lib64/perl5/5.10.0/x86_64-linux-thread-multi/CORE:/afs/desy.de/user/k/kiesel/scratch/singlePhoton/TreeWriter/
#export PATH=$PATH:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0/external/slc6_amd64_gcc491/bin/


export SCRAM_ARCH="slc6_amd64_gcc491"

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0
eval `scramv1 runtime -sh`
cd /afs/desy.de/user/e/eicker/BachelorTreeWriter
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
#LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0/biglib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0/lib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/llvm/3.6/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib:/opt/d-cache/dcap/lib64:/afs/desy.de/products/root/amd64_rhel60/5.34.00/lib:/usr/lib64/perl5/5.10.0/x86_64-linux-thread-multi/CORE:/afs/desy.de/user/k/kiesel/scratch/singlePhoton/TreeWriterNiklas
#PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0/bin/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0/external/slc6_amd64_gcc491/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/llvm/3.6/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/bin:/afs/desy.de/products/root/amd64_rhel60/5.34.00/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/bin:/usr/sge/bin:/usr/sge/bin/lx-amd64:/usr/lib64/qt-3.3/bin:/bin:/afs/desy.de/common/passwd:/usr/local/bin:/usr/bin:/afs/desy.de/user/k/kiesel/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin
#echo `pwd`
#ls
#echo $LD_LIBRARY_PATH
./treeWriter $*
echo Job ended sucessfully
