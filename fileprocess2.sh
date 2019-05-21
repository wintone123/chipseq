path=/mnt/c/chipseq/test8
dir_list="h3k4me1 h3k4me3 h3k9me3 h3k36me3"
chr_list="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"
cd $path

# for dir in $dir_list
# do
#     echo "spliting "$dir" ......"
#     for chr in $chr_list
#     do
#         echo "making chr "$chr" ......"
#         awk -va=$chr '{if($1 == a && $5 >= 40) print $0}' $dir"/"$dir".bed" > $dir"/"$dir"_chr"$chr".bed"
#     done
# done

for dir in $dir_list
do
    echo "spliting "$dir" ......"
    for chr in $chr_list
    do
        echo "making chr "$chr" ......"
        cat $dir"/"$dir".bed" | parallel --pipe awk -va=$chr \'{if\(\$1 == a \&\& \$5 \>= 40\) print \$0}\' | awk -va=$chr '{if($1 == a && $5 >= 40) print $0}' > $dir"/"$dir"_chr"$chr".bed"
    done
done