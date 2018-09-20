# Short-read WGS

Canonical SVs were called from WGS data. Details:

| Structural Variant | Chrom | Start    | End      | Variant caller |
| ------------------ |:-----:| :-------:| :------: | --------------:|
| Duplication 1      | X     | 17792578 | 18073738 | Canvas         |
| Inversion          | X     | 18074005 | 18532312 | Manta          |
| Duplication 2      | X     | 18248826 | 18532221 | Canvas         |

First, we will visualise the short-read sequencing data for the child and the two unaffected parents. Open IGV:

```
igv&
```

And open the following bam files (**Hint**: File > Load from file):

```
/mnt/albasj/data/illumina/child.chrX.bam
/mnt/albasj/data/illumina/father.chrX.bam
/mnt/albasj/data/illumina/mother.chrX.bam
```

Go to region:
chrX:17,699,000-18,653,000	

If you have problems visualising this region you can change the visibility range threshold to 1000 Kb in View>Preferences>Alignments.
If you still have problems, go to individual breakpoints and zoom out.

- Which is the mechanism of inheritance of this variant?
- Can you phase the cxSV from the short reads? Some informative SNPs are chrX:18064030, chrX:18073775, chrX:18504017.
