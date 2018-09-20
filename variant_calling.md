# Variant calling

Variants are called and stored in [VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf) format. This contains a header, and then data lines each containing information about a position in the genome.

<img src="//raw.githubusercontent.com/alsanju/train_malta_nanopore/master/images/vcf.png" alt="img_3" class="inline"/>

Currently, there are different algorithms for calling SVs from long-read sequencing data, including:
-	[Sniffles](http://github.com/fritzsedlazeck/Sniffles): best used with NGMLR. 
-	[NanoSV](http://github.com/philres/ngmlr): best used with LAST.

Since we used NGMLR for the alignment, now we will use sniffles for calling structural variants.

```
sniffles -m alignment/child.nanopore.ROI.sort.bam -v variant_calling/child.nanopore.ROI.vcf
```

How many SVs have been called?:

```
bcftools view -H variant_calling/child.nanopore.ROI.vcf | wc -l
```

The -s parameter can be changed to 1, where s is the minimum number of reads that support a SV (by default is 10).

```
sniffles -m alignment/child.nanopore.ROI.sort.bam -v variant_calling/child.nanopore.ROI.vcf -s 1
```

The information that is provided in sniffles’s output can be found in:
[http://github.com/fritzsedlazeck/Sniffles/wiki/Output](http://github.com/fritzsedlazeck/Sniffles/wiki/Output)


Finally, you can convert the VCF to a tab format:

```
scripts/vcf2tab.py variant_calling/child.nanopore.ROI.vcf
```

and inspect the SVs in IGV.
