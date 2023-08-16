module load samtools

indir=$1

for i in $1/*.vcf.gz; do
    a=$(bcftools query -f"%REF%ALT\n" $i | grep "^G\|^C" | wc -l)
    b=$(bcftools query -f"%REF%ALT\n" $i | wc -l)
    c=$(bcftools query -f"%REF%ALT\n" $i | grep "G$\|C$" | wc -l)
    #file_name GC_cont_REF GC_cont_ALT
    echo "$i $a $b $c" | awk '{print $2/$3 "\t" $4/$3}'
done > GC_prop
