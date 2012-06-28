indices=$2
names=( `cat ../data/shrtnames.txt`)
bisect=( `cat ../data/bisect_fine.txt`)
S=( `cat ../data/S_vals.txt`) 
N=( `cat ../data/N_vals.txt`) 

sadType="meteSAD empirSAD"

for i in $indices
do
  for j in $sadType
  do
    if [ $1 == unc ]
    then
       if [ $j == meteSAD ]
       then
         bsub -q week -J $i -o ./log_files/error_sim_${names[$i]}_$j.log\
           python spat_community_generation.py ${S[$i]} $N[$i]} 200 ${bisect[$i]}\
           False None ${names[$i]} 
       fi
       if [ $j == empirSAD ]
       then
         bsub -q week -J $i -o ./log_files/error_sim_${names[$i]}_$j.log\
           python spat_community_generation.py ${S[$i]} $N[$i]} 200 ${bisect[$i]}\ 
           False ../data/${names[$i]}_sad.csv ${names[$i]}
       fi
    fi
    if [ $1 == usu ]
    then
      if [ $j == meteSAD ] 
      then
        python spat_community_generation.py ${S[$i]} ${N[$i]} 10 ${bisect[$i]}\
        False None ${names[$i]}\
        > log_files/error_sim_${names[$i]}_$j.log 2>&1 &
      fi
      if [ $j == empirSAD ]
      then
        python spat_community_generation.py ${S[$i]} ${N[$i]} 200 ${bisect[$i]}\ 
        False ../data/${names[$i]}_sad.csv ${names[$i]}\
        > log_files/error_sim_${names[$i]}_$j.log 2>&1 &            
      fi
    fi
  done
done
