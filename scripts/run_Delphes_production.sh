#!/bin/bash

chmod 777 *.*

#set -o verbose

echo "Starting job on " `date`
echo "Running on " `uname -a`
echo "System release " `cat /etc/redhat-release`

export PRODHOME=`pwd`
export SCRAM_ARCH=slc6_amd64_gcc472
export RELEASE=CMSSW_5_3_19
export PATH=`pwd`:${PATH}


scram project ${RELEASE}
cd ${RELEASE}/src
eval `scram runtime -sh`
pwd

ln -s `which gcc` gcc
ln -s `which gfortran` f77
ln -s `which gfortran` g77

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
tar -xzf $PRODHOME/BRIDGE_MOD.tar.gz

tar -xzf $PRODHOME/Delphes.tar.gz
cd Delphes
chmod -R 777 *
#make clean
./configure
make
cd ..
tar -xzvf $PRODHOME/Cards.tar.gz
chmod -R 777 Cards/
cp Cards/MinBias_TuneZ2star_14TeV_pythia6_cff.py . 
cp BRIDGE_MOD/ApplyKfactor.py .
cp BRIDGE_MOD/kfactors_14TeV.dat .
cp BRIDGE_MOD/kfactors.dat .
cp BRIDGE_MOD/*.exe .
#wget http://red-gridftp11.unl.edu/Snowmass/HTBinned/CMSTP/14TEV/tt-4p-1100-1700-v1510_14TEV/tt-4p-1100-1700-v1510_14TEV_133094510.lhe.gz
cp -r BRIDGE_MOD/input .

# And then when all is ok

originDir=`pwd`
eosOutDir=/eos/cms/store/group/upgrade/delphes/alexis_tests/
#eosOutDir=/eos/cms/store/group/upgrade/delphes/ProdJun14/
export localDir=$PWD

startTime=`date +%s`

#input=tt-4p-1100-1700-v1510_14TEV_133094510.lhe.gz
theDir=B-4p-0-1-v1510_14TEV
input=B-4p-0-1-v1510_14TEV_102890512.lhe.gz
kFactor=kfactors_14TeV.dat
detCard=JetStudies_Phase_II_140PileUp_conf4.tcl
workDir=0
isBackground=1



if [ "$GLIDEIN_ResourceName" == USCMS-FNAL-WC1 ]; then
    export OSG_WN_TMP=${OSG_WN_TMP:-$PWD} # make sure that this is set
    cd $OSG_WN_TMP
fi

if [ -d ${workDir} ]; then rm -rf ${workDir}; fi
mkdir ${workDir}
cd ${workDir}


MadgraphDir=MG5
DelphesDir=Delphes
pythiaCard=pythia_card.dat

energy=14
DelphesVersion=`basename $DelphesDir`
nPU=`echo $detCard | cut -d '_' -f 4 | cut -d '.' -f 1`
process=`echo $input | cut -d '_' -f 1-2`
#seed=`echo $input | cut -d '_' -f 3 | cut -d '.' -f 1`
seed=102890517
part=5
total=5
phase=`echo $detCard | cut -d '_' -f 3`
if [ $phase == I ]; then
    configuration=0
else
    #configuration=`echo $detCard | cut -d '_' -f 5 | cut -d '.' -f 1 | cut -d 'f' -f 2`v2
    configuration=`echo $detCard | cut -d '_' -f 5 | cut -d '.' -f 1 | cut -d 'f' -f 2`
fi

#if [[ $isBackground == True ]]; then
#    DelphesOutput=${process}_${seed}_Phase${phase}_Conf${configuration}_${nPU}_seed${seed}_$part.root
#else
    DelphesOutput=`echo $input | cut -d '.' -f 1`_Phase${phase}_Conf${configuration}_${nPU}_seed${seed}_${part}of${total}.root
#fi
metaData=`echo $DelphesOutput | sed s/root/txt/`



cd ${localDir}/${workDir}
mkdir run
runDir=$PWD/run

mkdir Cards

#submitDir=$PWD
#MadgraphDir=$PWD/$MadgraphDir
MadgraphDir=$originDir/$MadgraphDir
#DelphesDir=$PWD/$DelphesDir
DelphesDir=$originDir/$DelphesDir
#BridgeDir=$PWD/${bridge}
BridgeDir=$originDir/${bridge}


cd $runDir

setupTime=`date +%s`

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Run LHE file through Bridge and Pythia for fragmentation + hadronization

cd $BridgeDir
echo $input
#i=`expr index "$input" _`
#theDir=${input:0:$i}14TEV

echo trying to copy file......
ls  ${runDir}

echo xrdcp root://eoscms.cern.ch//eos/cms/store/group/upgrade/lhe/${theDir}/${input} ${runDir}
xrdcp root://eoscms.cern.ch//eos/cms/store/group/upgrade/lhe/${theDir}/${input} ${runDir}

lhefile=${runDir}/unweighted_events.lhe
gunzip -c ${runDir}/$input > $lhefile
nEventsIn=`grep -c '<event>' $lhefile`
ls $lhefile
echo ==================================================== There are $nEventsIn for $lhefile ===============================================================

#qCut=`zgrep ' = xqcut' $originDir/$input  | awk '{print $1}'`
qCut=`zgrep ' = xqcut' ${runDir}/$input  | awk '{print $1}'`
if [ $qCut -ge "81" ]
    then
    exit 12
fi




echo "      QCUT = $qCut" >> ${localDir}/${workDir}/Cards/${pythiaCard}
if [ $process == "B-4p-0-1-v1510_14TEV" ]
    then
    echo >> ${localDir}/${workDir}/Cards/${pythiaCard}
    echo "IEXCFILE = 1" >> ${localDir}/${workDir}/Cards/${pythiaCard}
fi

echo $runDir
chmod 777 $lhefile
cp $lhefile $runDir
cd $runDir

maxjets=`grep "maxjetflav" $lhefile | awk '{print $1}'`
cp $PRODHOME/mgPostProcv2.py .
echo python ./mgPostProcv2.py -o file_out.lhe  -j ${maxjets} -q ${qCut} -e 5 -s $lhefile
python ./mgPostProcv2.py -o file_out.lhe  -j ${maxjets} -q ${qCut} -e 5 -s $lhefile

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
####### Apply kFactor is needed

#python $runDir/ApplyKfactor.py $runDir/$kFactor $lhefile ${runDir}/output_kfactors.lhe
#echo "Check for errors here,  sometimes the kfactor files has issues"
#$runDir/runDGE.exe n ${runDir}/input/particles.dat ${runDir}/input/SMparams.txt ${runDir}/input/interactions.dat ${runDir}/output_kfactors.lhe ${runDir}/unweighted_events.lhe ./grid/ time 2 finalstateparticles.txt

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
####### Run pythia via CMSSW with the "8TeV" hadronizer which includes tauola as well


cp $PRODHOME/hadronizer.py hadronizer.py

cmsRun hadronizer.py > pythia.log 2>&1

#if [ $? -ne 0 ]; then
#    cat pythia.log
#    cat $localDir/${workDir}/Cards/$pythiaCard
#    exit 14
#fi
cat pythia.log
#Pythiaxsec=`python $originDir/Cards/GetPythiaXS.py pythia.log  | grep sigma_matched  | awk '{print $3}'`
#new pythia xsec for CMSSW running
Pythiaxsec=` cat pythia.log  | grep "All included sub"  | awk '{print $10}'`
unit=` cat pythia.log | grep "\-\-\-\-\-\-\-\-\-\-\-\-\-\-I"  | awk '{print $2}'`

echo !!!!!!!!!!!!!!!!  We have the pythia step doneee !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! In $PWD
ls
chmod 777 -R $runDir

echo ===========================================================================
echo ===========================================================================
echo ===========================================================================
echo ===========================================================================
echo - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 0 - 


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#################Produce MinBias File

if [ $nPU != "NoPileUp" ]; then
    sed -i 's/50000/50/g' $originDir/Cards/MinBias_TuneZ2star_14TeV_pythia6_cff.py
    cmsRun  $originDir/Cards/MinBias_TuneZ2star_14TeV_pythia6_cff.py
    mv GenEvent_ASCII.dat GenEvent_ASCII.hepmc

    $DelphesDir/hepmc2pileup MinBias.pileup GenEvent_ASCII.hepmc
    if [ $? -ne 0 ]; then
	exit 15
    fi

    rm -rf GenEvent_ASCII.hepmc
fi

MBTime=`date +%s`


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#run Pythia output through Delphes

#ln -s ${runDir}/MinBias.pileup .
#$DelphesDir/DelphesSTDHEP $originDir/Cards/$detCard $DelphesOutput pythia_events.hep
sed -i "s/SET_SEEDHERE/$seed/g" $originDir/Cards/$detCard

$DelphesDir/DelphesCMSFWLite $originDir/Cards/$detCard $DelphesOutput file_out.root
if [ $? -ne 0 ]; then
    exit 16
fi

nEventsOut=`root -l -b -q $originDir/getNEvents.C\(\"$DelphesOutput\"\) | tail -n 1 | cut -d ')' -f 2` 

DelphesTime=`date +%s`

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#metadata

echo "User: " `eval whoami` > $metaData
echo "Date: " `date` >> $metaData
echo >> $metaData

echo "Process: " $process >> $metaData
echo "Pileup Conditions: " $nPU >> $metaData
echo "Phase: " $phase >> $metaData
echo "Configuration: " $configuration >> $metaData
echo "Random Seed: " $seed >> $metaData
echo >> $metaData

echo "Input LHE: " $input >> $metaData
echo "Input Events: " $nEventsIn >> $metaData
echo >> $metaData

echo "Pythia xsec in " $unit ":" $Pythiaxsec >> $metaData
echo "QCUT: " $qCut >> $metaData
echo >> $metaData

echo "Delphes Output: " $DelphesOutput >> $metaData
echo "Output Events: " $nEventsOut >> $metaData
echo "Delphes Version: " $DelphesVersion >> $metaData
echo "Detector Card: " $detCard >> $metaData
echo >> $metaData

echo "Minutes spent setting up job: " `expr $setupTime / 60 - $startTime / 60` >> $metaData
echo "Minutes spent running Pythia: " `expr $PythiaTime / 60 - $setupTime / 60` >> $metaData
echo "Minutes spent generating minBias events: " `expr $MBTime / 60 - $PythiaTime / 60` >> $metaData
echo "Minutes spent running Delphes: " `expr $DelphesTime / 60 - $MBTime / 60` >> $metaData
echo >> $metaData

echo "Pythia Card:" >> $metaData
cat ${localDir}/${workDir}/Cards/${pythiaCard} >> $metaData
echo >> $metaData

echo "Delphes Card:" >> $metaData
cat ${originDir}/Cards/$detCard >> $metaData
echo >> $metaData

echo "Condor Executable:" >> $metaData
cat $originDir/condor_exec.exe >> $metaData

cat $metaData

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Copy output

#change to eos asap
#cp ${DelphesOutput} $originDir
#cp ${metaData} $originDir
#copy=0


ls ${metaData} >> files_stageout

#command="xrdcp "${metaData}" root://eoscms.cern.ch/"${eosOutDir}""     


ls ${DelphesOutput} >> files_stageout

#cp $PRODHOME/copy.sh .
#chmod 777 copy.sh
#. copy.sh files_stageout $theDir

copy=0

#if [ $copy -ne 0 ]; then
#    exit 17
#fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
endTime=`date +%s`
echo "Time spent copying output (s): " `expr $endTime - $DelphesTime`
echo "Total runtime (m): " `expr $endTime / 60 - $startTime / 60`

#exit 0
