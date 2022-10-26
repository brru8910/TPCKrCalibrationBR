#!/bin/bash

if [[ $# -ne 2 ]]
then
    echo "Incorrect usage! Usage: ./runKryptonAnalyzerOnCondor.sh [Output Prefix] [Directory Containing Krypton.root Files]"
    exit 1
fi

exeName=`pwd`'/KryptonAnalyzer'

#Transfer arguments.
outputPrefix=$1
inputDirectory=$2

kryptonFileMatchPattern='krCalibration.root'

#Set up SHINE.
source /afs/cern.ch/user/n/na61qa/SHINEInstallations/runScript.sh

echo 'Running Krypton analysis on HTCondor.'
echo '[INFO] Output prefix: '$outputPrefix
echo '[INFO] Input directory: '$inputDirectory
echo '[INFO] Krypton file match pattern: '$kryptonFileMatchPattern
echo '[INFO] Condor log directory: '$kryptonFileMatchPattern

#De-activate glob.
#set -o noglob

arguments="-o "$outputPrefix" -i `ls "$inputDirectory"/*"$kryptonFileMatchPattern"`"
echo 'Exe command: '$exeName
# $arguments

#Create log directory.
timestamp=$(date +%Y%m%d_%H%M%S)
condorLogDirectory=`pwd`'/condorLogs/'$timestamp
mkdir -p $condorLogDirectory

# espresso     = 20 minutes
# microcentury = 1 hour
# longlunch    = 2 hours
# workday      = 8 hours
# tomorrow     = 1 day
# testmatch    = 3 days
# nextweek     = 1 week

queue='testmatch'

#Create sub file.
echo 'executable 	        =       '$exeName > condor.sub
echo 'arguments 	        =       '$arguments >> condor.sub
echo '+JobFlavour	        =	'$queue >> condor.sub
echo 'log		        =	'$condorLogDirectory'/condor.log' >> condor.sub
echo 'output		        =	'$condorLogDirectory'/stdout.log' >> condor.sub
echo 'error	        	=	'$condorLogDirectory'/error.log' >> condor.sub
echo 'transfer_input_files	=	bootstrap.xml,Config.txt' >> condor.sub
echo 'on_exit_remove            =       (ExitBySignal == False) && (ExitCode == 0)' >> condor.sub
echo 'max_retries               =       1' >> condor.sub
echo 'requirements              =       Machine =!= LastRemoteHost' >> condor.sub
echo 'getenv                    =       True' >> condor.sub
echo 'queue' >> condor.sub

#Submit job.
condor_submit condor.sub

#Re-activate glob.
#set +o noglob

