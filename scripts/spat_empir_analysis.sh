npar=$2
nperm=$3
commName=$4

metricsToCalc="sorensen"

dataType="abu binary"

if [ $1 == unc ] 
then
  for i in $commName           
  do
    for j in $metricsToCalc 
    do 
      for k in $dataType 
      do
        if [ $k == abu ]
        then
          export OMP_NUM_THREADS=$npar
          bsub -q week -M 20 -n $npar -R "span[hosts=1]" -J $i\
          -o ./log_files/error_$j_$k_$i.log\
          Rscript spat_empir_analysis.R $i $j $k $nperm $npar
        else
          export OMP_NUM_THREADS=1
          bsub -q week -M 20 -n 1 -J $i -o ./log_files/error_$j_$k_$i.log\
          Rscript spat_empir_analysis.R $i $j $k NULL 1
        fi
      done
    done
  done
fi

if [ $1 == usu ] 
then
  for i in $commName           
  do
    for j in $metricsToCalc 
    do 
      for k in $dataType 
      do
        if [ $k == abu ]
        then
          Rscript spat_empir_analysis.R $i $j $k $nperm $npar\
           >./log_files/error_$j_$k_$i.log 2>&1 &
        else
          Rscript spat_empir_analysis.R $i $j $k NULL 1\
          >./log_files/error_$j_$k_$i.log 2>&1 &
        fi
      done
    done
  done
fi
