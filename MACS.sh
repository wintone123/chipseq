path="/mnt/c/chipseq/test2"
cd $path

# MACS
for file in *_best.bed
do 
    basename=${file:0:-4}
    control="control_best.bed"
    if [ ${file:0:7} != 'control' ]
    then
        echo "calling peaks on "$file 
        macs -t $file -c $control -n $basename -f BED -g mm --nomodel
    fi
    echo "=============finish calling============="
done

# split file
chrom_list="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"
mkdir peaks_split_files
for file in *_peaks.bed
do
    basename=${file:0:-4}
    echo "spliting "$file
    for chr in $chrom_list
    do
        awk -va=$chr '{if($1 == "chr"a) print $0}' $file > "peaks_split_files/"$basename"_chr"$chr".bed"
    done
    echo "=============Done!============="
done
