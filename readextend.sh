# read extend
# Rscript readextend.R

# change file name   
path="/mnt/c/chipseq/test3"
cd $path

for file in ext_split_files/*_ext.txt
do 
    basename=${file:0:-4}
    echo "making "$basename".bed......"
    awk 'BEGIN{p=1}{if(p>1) print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$5; p++}' $file > $basename".bed"
done
rm ext_split_files/*_ext.txt
echo "==========done=========="