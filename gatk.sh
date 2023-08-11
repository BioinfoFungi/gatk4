mkdir -p output/index
mkdir -p output/mapping
mkdir -p output/base_recalibrator
mkdir -p output/vcf
mkdir -p output/gvcf

WORK_DIR=/home/wangyang/workspace/GATK/
ln -s $WORK_DIR/2-germline/ref/ref.fasta output/index



########################################################################################
# 比对
########################################################################################
bwa index output/index/ref.fasta


bwa mem -t 10 -R '@RG\tID:foo_lane\tPL:illumina\tLB:library\tSM:mother' output/index/ref.fasta \
    2-germline/fastq/mother_1.fastq \
    2-germline/fastq/mother_2.fastq > output/mapping/mother.sam


samtools view -bS  -@ 10 output/mapping/mother.sam -o  output/mapping/mother.bam

# samtools  flagstat -@ 10 output/mapping/mother.bam

samtools  sort output/mapping/mother.bam -@ 10 -o   output/mapping/mother.sorted.bam


########################################################################################
# 去除PCR重复
########################################################################################


gatk MarkDuplicates -I output/mapping/mother.sorted.bam \
    -O output/mapping/mother.sorted.markdup.bam  \
    -M output/mapping/mother.markdup.metrics.txt


samtools index output/mapping/mother.sorted.markdup.bam


# samtools  flagstat -@ 10 output/mapping/mother.sorted.markdup.bam 
# samtools view -f 1024  output/mapping/mother.sorted.markdup.bam | wc -l
# samtools view -c output/mapping/mother.sorted.markdup.bam


########################################################################################
# 碱基质量矫正
########################################################################################
# 这里为了测试BQSR，直接使用测试数据提供的bam生成vcf文件作为群体中已知的indel，计算出了所有需要进行重校
# 正的read和特征值，然后把这些信息输出为一份校准表文件

gatk HaplotypeCaller \
    -R 2-germline/ref/ref.fasta \
    -I 2-germline/bams/mother.bam \
    -O output/vcf/motherHC.vcf \
    -L 20:10,000,000-10,200,000
bcftools  view output/vcf/motherHC.vcf | less -S




# bcftools view -H vcf-resource/Homo_sapiens_assembly38.known_indels.vcf.gz | cut -f 1 | sort | uniq
samtools faidx output/index/ref.fasta
gatk CreateSequenceDictionary -R output/index/ref.fasta -O output/index/ref.dict 


# bcftools view -h vcf-resource/Homo_sapiens_assembly38.known_indels.vcf.gz | grep "^##"  >  output/vcf/Homo_sapiens_assembly38.known_indels.chr22.vcf
# bcftools view -h vcf-resource/Homo_sapiens_assembly38.known_indels.vcf.gz | grep "#CHROM"  >>output/vcf/Homo_sapiens_assembly38.known_indels.chr22.vcf
# bcftools  view vcf-resource/Homo_sapiens_assembly38.known_indels.vcf.gz chr22   >> output/vcf/Homo_sapiens_assembly38.known_indels.chr22.vcf
# bgzip output/vcf/Homo_sapiens_assembly38.known_indels.chr22.vcf
# bcftools  view vcf-resource/Homo_sapiens_assembly38.known_indels.vcf.gz chr22 | bgzip > output/vcf/Homo_sapiens_assembly38.known_indels.chr22.vcf.gz
# tabix -p vcf output/vcf/Homo_sapiens_assembly38.known_indels.chr22.vcf.gz

gatk BaseRecalibrator \
    -I output/mapping/mother.sorted.markdup.bam  \
    -R output/index/ref.fasta \
    --known-sites output/vcf/motherHC.vcf \
    -O output/base_recalibrator/mother.recal_data.table

gatk ApplyBQSR \
    -R output/index/ref.fasta \
    -I output/mapping/mother.sorted.markdup.bam  \
    --bqsr-recal-file output/base_recalibrator/mother.recal_data.table \
    -O output/mapping/mother.sorted.markdup.BQSR.bam

samtools index output/mapping/mother.sorted.markdup.BQSR.bam

########################################################################################
# 变异检测-单样本模式
########################################################################################

gatk HaplotypeCaller \
    -R output/index/ref.fasta \
    -I output/mapping/mother.sorted.markdup.BQSR.bam \
    -O output/vcf/mother.sorted.markdup.vcf \
    -L 20:10,000,000-10,200,000
 

gatk HaplotypeCaller \
    -R output/index/ref.fasta \
    -I output/mapping/mother.sorted.markdup.BQSR.bam \
    -O output/vcf/mother.sorted.markdup.vcf \
    -D 2-germline/resources/dbsnp.vcf \
    -L 20:10,000,000-10,200,000
 


bcftools  view 2-germline/resources/dbsnp.vcf | less -S
bcftools  view output/vcf/mother.sorted.markdup.vcf  | less -S
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



vep_install -a cf -s homo_sapiens -y GRCh38 -c db/vep --CONVERT


########################################################################################
# 变异检测-变异的质控
########################################################################################
gatk HaplotypeCaller \
    -R output/index/ref.fasta \
    -I output/mapping/mother.sorted.markdup.BQSR.bam \
    -O output/vcf/mother.sorted.markdup.vcf \
    -D 2-germline/resources/dbsnp.vcf \
    -L 20:10,000,000-10,200,000
 
# 使用SelectVariants，选出SNP
gatk SelectVariants \
    -select-type SNP \
    -V output/vcf/mother.sorted.markdup.vcf  \
    -O output/vcf/mother.sorted.markdup.snp.vcf 

gatk VariantFiltration \
    -V output/vcf/mother.sorted.markdup.snp.vcf  \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O output/vcf/mother.sorted.markdup.snp.filter.vcf 


# 使用SelectVariants，选出Indel
gatk SelectVariants \
    -select-type INDEL \
    -V output/vcf/mother.sorted.markdup.vcf \
    -O output/vcf/mother.sorted.markdup.indel.vcf 

gatk VariantFiltration \
    -V output/vcf/mother.sorted.markdup.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O ../output/E.coli/E_coli_K12.indel.filter.vcf.gz





