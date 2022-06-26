#!/usr/bin/env bash
##File for output xs files on demand by flags

##read flags (can be only one letter)
while getopts a:b:c:d:e:f:g: flag
do
    case "${flag}" in
        a) theta1=${OPTARG};;
        b) phi1=${OPTARG};;
        c) theta2=${OPTARG};;
	d) phi2=${OPTARG};;
	e) interaction=${OPTARG};;
	f) ann=${OPTARG};;
	g) energy=${OPTARG};;
    esac
done

##energy="13MeV"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
indirectory="${SCRIPT_DIR}/scattering_lengths/${energy}/${ann}/${interaction}"
outdirectory="${indirectory}/theoretical_xs" #"${SCRIPT_DIR}/${interaction}/theoretical_xs"
infile="${indirectory}/input1.dat"
fortran_file="CDBD4u1E13.0_3N"


##create new subdirectoy in execution directory if it does not exist
mkdir -p ${outdirectory}

outfile="${outdirectory}/xs_${theta1}_${phi1}_${theta2}_${phi2}.dat"

##echo "xs for: ${theta1}    ${phi1}   ${theta2}    ${phi2}"
##gawk -i inplace -v tt1="${theta1}" -v tt2="${theta2}" -v pt1="${phi1}" -v pt2="${phi2}" '(NR==12){$1=tt1; $2=pt1; $3=tt2; $4=pt2}1' $infile
##conservative way of doing awk since IGFAE's server does not have gawk -i
awk -v tt1="${theta1}" -v tt2="${theta2}" -v pt1="${phi1}" -v pt2="${phi2}" '(NR==12){$1=tt1; $2=pt1; $3=tt2; $4=pt2}1' $infile > tmp && mv tmp ${infile}
##cd because BreakUpObs needs to be run in that directory
cd "${indirectory}"
./BreakUpObs.x-cs < $fortran_file > $outfile
##this does not do anything! But as long as we use global directory names, it does not matter
cd ${SCRIPT_DIR}
##and now use sed to delete first 4 lines in each file (we also can have used awk )
sed -i '1,4d' $outfile
