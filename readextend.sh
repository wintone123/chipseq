# read extend
# Rscript readextend.R

# change file name   
path="/mnt/c/chipseq/test3"
cd $path

for file in *.txt
do 
    basename=${file:0:-4}
    awk 'BEGIN{p=1}{if(p>1) print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$5; p++}' $file > $basename"_ext.bed"
done
echo "==========done=========="