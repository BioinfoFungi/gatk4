
bcftools mpileup -f \
    2-germline/ref/ref.fasta \
    2-germline/bams/mother.bam  \
    | bcftools call -mv -o output/bcftools/mother.call.vcf 


bgzip output/bcftools/mother.call.vcf 
#默认是构建.csi索引
bcftools index output/bcftools/mother.call.vcf.gz
# 另外还有一种.tbi索引
bcftools index -t output/bcftools/mother.call.vcf.gz





