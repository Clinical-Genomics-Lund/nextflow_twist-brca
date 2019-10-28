# Sentieon QC script version 1.01
git clone https://github.com/Clinical-Genomics-Lund/qc_sentieon.git tmp
git --git-dir=./tmp/.git --work-tree=./tmp checkout v1.02
cp tmp/qc_sentieon.pl bin/.
rm -rf tmp/
# Aggregate VCF script version 1.01
git clone https://github.com/Clinical-Genomics-Lund/aggregate_vcf.git tmp
git --git-dir=./tmp/.git --work-tree=./tmp checkout v1.02
cp tmp/aggregate_vcf.pl bin/.
rm -rf tmp/
