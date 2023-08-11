mkdir -p output/index
mkdir -p output/mapping
mkdir -p output/base_recalibrator
mkdir -p output/vcf
mkdir -p output/gvcf

########################################################################################
# 变异检测-单样本模式
########################################################################################

gatk HaplotypeCaller \
    -R output/index/ref.fasta \
    -I 2-germline/bams/mother.bam \
    -I 2-germline/bams/father.bam \
    -O output/vcf/mutiple_sample.vcf \
    -L 0
 
bcftools  view output/vcf/mutiple_sample.vcf  | less -S
bcftools  view  output/vcf/mother.sorted.markdup.vcf  | less -S
########################################################################################
# 变异检测-多样本模式
########################################################################################
gatk HaplotypeCaller \
    -R output/index/ref.fasta \
    -I output/mapping/mother.sorted.markdup.BQSR.bam \
    -O output/gvcf/mother.sorted.markdup.g.vcf \
    --emit-ref-confidence GVCF \
    -L 20:10,000,000-10,200,000


gatk CombineGVCFs \
    -R output/index/ref.fasta  \
    --variant output/gvcf/mother.sorted.markdup.g.vcf  \
    -O output/gvcf/cohort.g.vcf.gz
 
gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R output/index/ref.fasta \
    -V output/gvcf/cohort.g.vcf.gz \
    -O output/gvcf/output.vcf.gz
 
bcftools  view output/gvcf/output.vcf.gz | less -S