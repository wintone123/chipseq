# remove low quality reads
for file in *.bed
do
    basename=${file:0:-4}
    echo "creating "$basename"_best.bed"
    awk '{if($5 >= 40) print $0}' $file > $basename"_best.bed" 
    echo "==========Done!=========="
done

# split file as chromosome
mkdir split_files
chrom_list="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"
for file in *_best.bed
do
    basename=${file:0:-4}
    echo "spliting "$file
    for i in $chrom_list
    do
        awk -va=$i '{if($1 == "chr"a) print $0}' $file > "split_files/"$basename"_chr"$i".bed"
    done
    echo "==========Done!=========="
done