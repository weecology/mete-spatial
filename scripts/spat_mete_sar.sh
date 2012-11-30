indices=$2
names=( `cat ../data/shrtnames.txt`)

for i in $indices
do
  if [ $1 == unc ]
    bsub -q week -M 8 -J ${names[$i]}\
    -o ./log_files/error_mete_sar_${names[$i]}.log\ 
    python spat_mete_sar.py ${names[$i]}\
  fi
  if [ $1 == usu ]
    python spat_mete_sar.py ${names[$i]}\
    > ./log_files/error_mete_sar_${names[$i]}.log 2>&1 &  
  fi
done
