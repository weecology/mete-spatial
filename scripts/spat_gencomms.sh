for S in 10  11  13  14  16  18  21  23  26  30  34  38  43  48  55  62  70  78  89 100



do
 
 
  for N in 120 186 289 447 694 1076 1668 2587 4011 6220 9646 14957 23193 35965 55769 86479 134099 207941 322444 500000


  do


    python spat_community_generation.py $S $N 200 13 False None space &

    sleep 20

  done
done


