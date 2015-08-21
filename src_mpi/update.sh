

for F in * 
do
done




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


