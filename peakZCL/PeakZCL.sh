#!/bin/bash

printf "\n\n===================================================================\n PeakZCL. A peak caller tool based on Zero-Crossing Lines for ChIP-seq\n===================================================================\n\n"


day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "Started at $day $hour \n"

sPath="`dirname \"$0\"`"
sPath="`( cd \"$sPath\" && pwd )`"

#printf $sPath\n

usage()
{
cat << EOF

OPTIONS:
   -f   Directory containing Sample file (required)
   -c   Directory containing control file
   -n   Normalization method ("SES","SDS", "N", default: SES)
   -w   Decimating signal value (default: 50)
   -l   Length signal wavelet (default: 1000000)
   -s   Scales of the wavelet (default: 30)
   -g   Zero-crossing lines leverage threshold (default: 5)
   -i   Baseline intensity threshold (default: 0)
   -j   Bp range between peaks to clustering purposes (default: 500)
   -a   Minimum value of the area of the selected peaks given in log10(x) (default:1.6)
   -k   Difference between peak area of sample and input. The value is given in log10(x) (default :0.2)
   -o   Output BED file name (default: ZCL_peaks.final.bed)
   -v   Path to annotation database in .gtf format
   -p   Number of processors used by R scripts (default: 1)
EOF
}

#Default parameters for broad peak signals --

normalization="SES"
decimating="50"
length="1000000"
scales="30"
levth="5"
threshold="0"
clustering="500"
area="1.6"
foldchange="0.2"
cores="1"
sdir=""
cdir=""
odir="NULL"
outnamefile="ZCL_peaks"
annotfile=""

while getopts "f:c:n:w:l:s:g:i:j:a:k:o:v:y:p:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
	;;
	c) cdir=$OPTARG
	;;
	n) normalization=$OPTARG
	;;
	w) decimating=$OPTARG
	;;
	l) length=$OPTARG
	;;
	s) scales=$OPTARG
	;;
	g) levth=$OPTARG
	;;
	i) threshold=$OPTARG
	;;
	j) clustering=$OPTARG
	;;
	a) area=$OPTARG
	;;
	k) foldchange=$OPTARG
	;;
	o) outnamefile=$OPTARG
	;;
	v) annotfile=$OPTARG
	;;
  p) cores=$OPTARG
  ;;
	?)
	usage
	exit
	;;
	esac
done

if [[ -z $sdir ]]
then
     usage
     exit 1
fi

printf "preprocessing chip and input files\n"


EXP=$sdir'chr_names.txt'


mkdir $sdir'TMP'
mkdir $sdir'TMP/chip/'

if [[ ! -z $cdir ]]; then
  mkdir $sdir'TMP/input/'

fi
#wdir= $sdir/TMP

#if [[ $odir == "NULL" ]]; then
mkdir $sdir'OUTPUT/'
odir=$sdir'OUTPUT/'

#fi
export LANG=C
export LC_ALL=C

# preprocessing chip
if [[ ! -z $sdir ]]; then
  printf "creating intermediate chip files\n"

  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
    chrom=$line
    printf "$line\n"
    grep "$chrom\>" $sdir/*.bedGraph >> $sdir/TMP/chip/$line.chrom &
  done < "${EXP}"

fi
wait
#preprocessing input

if [[ -d $sdir/TMP/input ]]; then
  printf "creating intermediate input files\n"
  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
    chrom=$line
    grep "$chrom\>" $cdir/*.bedGraph >> $sdir/TMP/input/$line.chrom &
  done < "${EXP}"
fi
wait

printf "DONE!\n"
##launch Rscript for each chromosome

#humanchrlist=( $EXP )
printf "$humanchrlist"\n

if [[ ! -z $cdir ]]; then
  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

               outnamefileR=$line
               printf "Running Rscript..... $sPath/ZCL_CHIP.R -f $sdir/TMP/chip/$line.chrom -c $sdir/TMP/input/$line.chrom -d "$odir" -n "$normalization" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -i "$threshold" -a "$area" -k "$foldchange" -o "$outnamefileR" \n" &
               Rscript $sPath/ZCL_CHIP.R -f $sdir/TMP/chip/$line.chrom -c $sdir/TMP/input/$line.chrom -d "$odir" -n "$normalization" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -i "$threshold" -a "$area" -k "$foldchange" -o "$outnamefileR" &
               NPROC=$(($NPROC+1))
               if [ "$NPROC" -ge "$cores" ]; then
                 wait
                 NPROC=0
               fi

  done < "${EXP}"
else
    while IFS='' read -r line || [[ -n "$line" ]]; do
              outnamefileR=$line
              printf "Running Rscript ..... $sPath/ZCL_CHIP.R -f $sdir/TMP/chip/$line.chrom -d "$odir" -n "$normalization" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -i "$threshold" -a "$area" -k "$foldchange" -o "$outnamefileR" \n" &
              Rscript $sPath/ZCL_CHIP.R -f $sdir/TMP/chip/$line.chrom -d "$odir" -n "$normalization" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -i "$threshold" -a "$area" -k "$foldchange" -o "$outnamefileR" &
              NPROC=$(($NPROC+1))
              if [ "$NPROC" -ge "$cores" ]; then
                wait
                NPROC=0
              fi
  done < "${EXP}"
fi

wait

printf "concatenate output files\n"
#concatenate files
cat $odir*.bed >> $odir$outnamefile.bed
#sort -k1,1V -n $odir$outnamefile.bed >> $odir$outnamefile.final.bed
#rm $odir$outnamefile.bed
Rscript $sPath/process_Output_bed.R $odir$outnamefile.bed

peaks= wc -l $odir$outnamefile"_ordered".bed
#rm $odir$outnamefile"_ordered".bed
printf "$peaks peaks found\n"

if [[ ! -z $annotfile ]]; then
  Rscript $sPath/Annotation.R $annotfile $odir$outnamefile"_ordered".bed $odir$outnamefile"_annot.txt"

fi




rm -rf $sdir/TMP

day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "finished at $day $hour !\n"
