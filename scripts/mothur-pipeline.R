# mothur pipeline

# done
make.file(inputdir=., type=gz, prefix=stability)
make.contigs(file=stability.files, processors=8)
summary.seqs(fasta=current, count=current)
screen.seqs(fasta=current, count=current, maxambig=0, maxlength=430, maxhomop=8)
summary.seqs(fasta=current, count=current)
unique.seqs(fasta=current, count=current)
summary.seqs(count=current)
align.seqs(fasta=current, reference=silva.nr_v132.v34.unique.align)
summary.seqs(fasta=current, count=current)
screen.seqs(fasta=current, count=current, start=2, end=17012)
summary.seqs(fasta=current, count=current)
filter.seqs(fasta=current, vertical=T, trump=.)
unique.seqs(fasta=current, count=current)
pre.cluster(fasta=current, count=current, diffs=2)
chimera.vsearch(fasta=current, count=current, dereplicate=t)
summary.seqs(fasta=current, count=current)
classify.seqs(fasta=current, count=current, reference=trainset18_062020.pds.fasta, taxonomy=trainset18_062020.pds.tax)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=current, count=current)
rename.file(fasta=current, count=current, taxonomy=current, prefix=final)
cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=4, cutoff=0.03)

# to do