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
        bsub -q week -J $commName\
        Rscript spat_empir_analysis.R $i $j $k\
        -o "error_${j}_${k}_${i}.log"
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
        Rscript spat_empir_analysis.R $i $j $k\
        >"error_${j}_${k}_${i}.log" 2>&1 &
      done
    done
  done
fi