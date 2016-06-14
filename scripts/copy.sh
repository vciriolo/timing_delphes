theDir=$2
while read line
do
echo will now copy  xrdcp $line root://eoscms.cern.ch//eos/cms/store/group/upgrade/delphes/ProdJun14/$theDir/$line
xrdcp $line root://eoscms.cern.ch//eos/cms/store/group/upgrade/delphes/ProdJun14/$theDir/$line
#xrdcp $line root://eoscms.cern.ch//eos/cms/store/group/upgrade/delphes/alexis_tests/DelphesJun20_v2/$theDir/$line
done<$1
