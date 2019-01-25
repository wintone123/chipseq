path="/mnt/c/chipseq/test3"
cd $path

# sam to bam
for file in *.sam
do
    echo "transfering "$file"......"
    basename=${file:0:-4}
    samtools view -bS -q 40 file > $basename".bam"
    echo $basename".bam established!"
done 

# bam to bad
for file in *.bam
do
    echo "transfering "$file"......"
    basename=${file:0:-4}
    bedtools bamtobed -i file > $basename".bed"
    echo $basename".bed established!"
done

haha