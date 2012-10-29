commName=$2

metricsToCalc="sorensen varWithin"

dataType="abu binary"

if [ $1 == unc ] 
then
  for i in $commName           
  do
    for j in $metricsToCalc 
    do 
      for k in $dataType 
      do
        bsub -q week -x -M 10 -J $i -o ./log_files/error_$j_$k_$i.log\
        Rscript spat_empir_analysis.R $i $j $k 
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
        Rscript spat_empir_analysis.R $i $j $k >./log_files/error_$j_$k_$i.log 2>&1 &
      done
    done
  done
fi
