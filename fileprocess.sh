path="/mnt/c/chipseq/test3"
cd $path

# chrom list
chrom_list="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"

# load environment
source activate bioinfo2

echo "=============Let's Go!============="

# sam to bam (remove quality < 40 reads)
# for file in *.sam
# do
#    echo "transfering "$file"......"
#    basename=${file:0:-4}
#    samtools view -bS -q 40 $file > $basename".bam"
#    echo $basename".bam established!"
# done 
# echo "=============Done!============="

# bam to bad
for file in *.bam
do
    echo "transfering "$file"......"
    basename=${file:0:-4}
    bedtools bamtobed -i $file > $basename".bed"
    echo $basename".bed established!"
done
echo "=============Done!============="

# remove low quality reads
for file in *.bed
do
    basename=${file:0:-4}
    echo "creating "$basename"_best.bed......"
    awk '{if($5 >= 40) print "chr"$0}' $file > $basename"_best.bed" 
done
echo "=============Done!============="

# bed file split
outdir="split_files"
mkdir $outdir
for file in *_best.bed
do
    basename=${file:0:-4}
    echo "spliting "$file"......"
    for i in $chrom_list
    do
        awk -va=$i '{if($1 == "chr"a) print $0}' $file > $outdir"/"$basename"_chr"$i".bed"
    done
done
echo "=============Done!============="

# peak calling by macs2
for file in *_best.bed
do 
    basename=${file:0:-4}
    control="control_best.bed"
    if [ $file != 'control_best.bed' ]
    then
        echo "calling peaks on "$file 
        macs2 callpeak -t $file -c $control -n $basename -f BED -g mm -B -q 0.01
    fi
done
echo "=============finish calling============="

# peak file split
outdir="narrowPeak_spilt_files"
mkdir $outdir
for file in *.narrowPeak
do
    basename=${file:0:-11}
    echo "spliting "$file
    for i in $chrom_list
    do
        awk -va=$i '{if($1 == "chr"a) print $0}' $file > $outdir"/"$basename"_chr"$i".narrowPeak"
    done
done
echo "=============Done!============="

# read extension
# Rscript readextend.R

# create new bed file
# for file in ext_split_files/*_ext.txt
# do 
#     basename=${file:0:-4}
#     echo "making "$basename".bed......"
#     awk 'BEGIN{p=1}{if(p>1) print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$5; p++}' $file > $basename".bed"
# done
# rm ext_split_files/*_ext.txt
# echo "=============Done!============="
