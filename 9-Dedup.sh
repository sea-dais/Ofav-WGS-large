idev -p spr -N 1 -n 1 -t 03:00:00 -A IBN21018

bcftools norm -d all vcf/7x11_ind3.vcf.gz -Oz -o vcf/7x11_ind3_dedup.vcf.gz
tabix -p vcf vcf/7x11_ind3_dedup.vcf.gz
bcftools view -H vcf/7x11_ind3_dedup.vcf.gz | cut -f1,2 | sort | uniq -d | wc -l   # -> 0

bcftools norm -d all vcf/11x7_ind3.vcf.gz -Oz -o vcf/11x7_ind3_dedup.vcf.gz
tabix -p vcf vcf/11x7_ind3_dedup.vcf.gz

# confirm the duplicates are gone (should now be 0)
bcftools view -H vcf/11x7_ind3_dedup.vcf.gz | cut -f1,2 | sort | uniq -d | wc -l