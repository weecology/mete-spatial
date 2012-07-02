indices=$2

for i in $indices
do
  if [ $1 == unc ]
  then
    bsub -q day -M 6 -J parm_space -o ./log_files/error_sim_S$i.log\
    python spat_gencomms.py $i
  fi
  if [ $1 == usu ] 
  then
    python spat_gencomms.py $i > ./log_files/error_sim_S$i.log 2>&1 &
  fi
done
