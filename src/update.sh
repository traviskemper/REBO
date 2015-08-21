

for F in * 
do
sed '/USE ANALYSIS/d' $F | sed '/USE BEAM/d' | sed '/USE MIN/d' > ../src-serialv3/$F
done

sed '/USE ANALYSIS/d' ../s$F | sed '/USE BEAM/d' | sed '/USE MIN/d' > ../src-serialv3/$F



 diff default.f varrm.txt | awk '{if ($1!='!') print $1}' 


rm  tochange.txt
for KW in `diff default.f varrm.txt | awk '{print $2}'`
do
 if [ "KW" != '!' ] ; then
  echo $KW >> tochange.txt
  grep -in "$KW" *.f  >>  tochange.txt
 fi 
done


for KW in `diff default.f varrm.txt | awk '{print $2}'`
do
 if [ "KW" != '!' ] ; then
  echo $KW >> tochange.txt
  grep -in "$KW" *.f  >>  tochange.txt
 fi
done

rm  tochange.txt
for KW in ` cat e.txt`
do
   echo $KW >>  tochange.txt 
   grep -in "$KW" *.f   >>  tochange.txt 
done


