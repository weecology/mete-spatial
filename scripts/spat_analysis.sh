indices=$2

## arguments for job
S=( `cat ../data/S_vals.txt`)
N=( `cat ../data/N_vals.txt`)
bisec_fine=( `cat ../data/bisect_fine.txt`)
bisec_coarse( `cat ../data/bisect_coarse.txt`)
grain_fine=( `cat ../data/grain_fine.txt`)
names=( `cat ../data/shrtnames.txt`)

sadType="meteSAD empirSAD"

dataType="abu binary"

metrics="sorensen varWithin"

for i in $indices
do
  for j in $sadType
  do
    for k in $dataType
    do
      for m in $metrics
      do
        if [ $1 == unc ]
        then
          if [ $j == meteSAD ]
          then
            bsub -q week -M 8 -J ${names[$i]}\
            -o ./log_files/error_sim_analysis_${names[$i]}_$j.log\
            Rscript spat_analysis.R ${S[$i]} ${N[$i]} 200\
            ${bisect_fine[$i]} ${bisect_coarse[$i]}\
            ${grain_fine[$i]} False $k $m NA NA ${names[$i]} TRUE 
          fi
          if [ $j == empirSAD ]
          then
            bsub -q week -M 8 -J ${names[$i]}\
            -o ./log_files/error_sim_analysis_${names[$i]}_$j.log\
            Rscript spat_analysis.R ${S[$i]} ${N[$i]} 200\
            ${bisect_fine[$i]} ${bisect_coarse[$i]}\
            ${grain_fine[$i]} False $k $m NA NA ${names[$i]}_$j TRUE  
          fi
        fi
        if [ $1 == usu ]
        then
          if [ $j == meteSAD ] 
          then
            Rscript spat_analysis.R ${S[$i]} ${N[$i]} 200\ 
            ${bisect_fine[$i]} ${bisect_coarse[$i]}\
            ${grain_fine[$i]} False $k $m NA NA ${names[$i]} TRUE\
            > log_files/error_sim_analysis_${names[$i]}_$j.log 2>&1 &
          fi
          if [ $j == empirSAD ]
          then
            Rscript spat_analysis.R ${S[$i]} ${N[$i]} 200\ 
            ${bisect_fine[$i]} ${bisect_coarse[$i]}\
            ${grain_fine[$i]} False $k $m NA NA ${names[$i]}_$j TRUE\
            > log_files/error_sim_analysis_${names[$i]}_$j.log 2>&1 &
          fi
        fi
      done
    done
  done
done
