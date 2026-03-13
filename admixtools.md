## admixtools

install from conda:
```
conda activate clarkii
conda install -c bioconda eigensoft
```

to install admixtools from source in my mac (helped by ChatGPT)

```
git clone https://github.com/DReichLab/AdmixTools.git
cd Admixtools/src
mv Makefile ..
mv Makefile.mac Makefile

nano compat_strchrnul.c
# paste this in file
#include <string.h>

char *strchrnul(const char *s, int c) {
    while (*s && *s != (char)c) {
        s++;
    }
    return (char *)s;
}
######

nano Makefile

# remove top of file and replace with this:
HOMEL=$(PWD)
TOP=../bin
BIN=$(HOMEL)/../bin

BREW := $(shell brew --prefix)

### *** needs argp. run brew install argp-standalone
override LDLIBS += compat_strchrnul.o $(BREW)/lib/libargp.a -lgsl -lblas -llapack -lm -lnick
# Some Linux distributions require separate lapacke library
# override LDLIBS += -llapacke

override LDFLAGS += -g -L./nicksrc -L$(BREW)/lib
override CFLAGS += -c -g -Wimplicit -I./ -I./nicksrc -I$(BREW)/include

# Harvard Medical School O2 cluster additions
ifdef SLURM_CONF
override CFLAGS += -I/n/app/openblas/0.2.19/include -I/n/app/gsl/2.3/include
override LDFLAGS += -L/n/app/openblas/0.2.19/lib -L/n/app/gsl/2.3/lib/
endif

ND = nicksrc
NLIB = $(ND)/libnick.a

#######

brew install argp openblas gsl argp-standalone
make clobber
make clean
cc -c -g compat_strchrnul.c -o compat_strchrnul.o
make
```

### prep input

filter out non "nc_" chromosomes (smaller ones):
```
bcftools view -h data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf.gz | grep -v "contig=<ID=NW_" > data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf
bcftools view -H data/vcfs/pclarkii.qc_ac_bial.lowmiss_maf.bcf.gz | awk '$1 ~ /^NC_/' >> data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf
```
make plink file
```
plink \
  --vcf data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf \
  --double-id \
  --allow-extra-chr \
  --chr-set 94 \
  --make-bed \
  --out data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr
```
make a chromosome to number map and change bim chromosome names to numbers:
```
# map
paste -d ' ' \
    <(grep -v "#" data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf | cut -f1 | uniq) \
    <(for n in {1..94}; do echo $n; done) \
    > data/admixtools/chrom_map.txt
# change names:
awk 'NR==FNR {map[$1]=$2; next} {$1=map[$1]; print}' data/admixtools/chrom_map.txt \
    data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.bim > data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bim
cp data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.bed data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bed
cp data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.fam data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.fam
```

from and Sonia's code to generate input:
```
nano data/admixtools/par.PLINK2EIGEN

# paste this in file:
genotypename:    data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bed
snpname:         data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.bim
indivname:       data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.num_chr.fam
outputformat:    EIGENSTRAT
genotypeoutname: data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.geno
snpoutname:      data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.snp
indoutname:      data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.ind
familynames:     NO
####### 

conda activate clarkii
convertf -p data/admixtools/par.PLINK2EIGEN
```

for ind file:
```
paste -d ' ' \
    <(bcftools query -l data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf) \
    <(yes "U" | head -n 328) \
    <(bcftools query -l data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.nc_chr.vcf | cut -c 1,2) \
    > data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.ind
```

for the snp file:
```
awk '{sub(/\./,"snp" NR)}1' data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.snp \
    > tmp && mv tmp data/admixtools/pclarkii.qc_ac_bial.lowmiss_maf.eigenstrat.snp
```
[admixr](https://cran.r-project.org/web/packages/admixr/vignettes/vignette-01-tutorial.html)