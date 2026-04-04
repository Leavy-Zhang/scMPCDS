#!/conda/bin/python

import os,re
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen




os.chdir('vcf_dir')

for i in os.listdir():
    os.system('bgzip %s' % i)
    os.system('tabix -p vcf %s.gz' % i)

os.system('bcftools merge *.vcf.gz -o merged.vcf.gz')

genInstall.install('GRCh38', rsync=False, bash=True)
matGen.SigProfilerMatrixGeneratorFunc("mitoAnalysis", "GRCh38", "/directory/of/merged_vcf_file", chrom_based=True, plot=True, tsb_stat=True )

