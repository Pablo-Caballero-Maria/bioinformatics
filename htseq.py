import HTSeq
import itertools
from matplotlib import pyplot

'''
bam_reader = HTSeq.BAM_Reader("SRR1039508_sorted.bam")
# this prints the first 5 reads in the bam file, (along with their sequence, quality, chromosome, and strand)
for alignment in itertools.islice(bam_reader, 5):
    print("Name: ", alignment.read.name)
    print("Sequence: ", alignment.read.seq)
    print("Quality: ", alignment.read.qual)
    print("Chromosome: ", alignment.iv.chrom)
    print("="*100)
'''

'''
coverage = HTSeq.GenomicArray("auto", stranded=True, typecode='i')
# this stores the coverage of the first 1000 reads in the bam file, the coverage increases if the read is aligned
i = 0
for alignment in itertools.islice(bam_reader, 1000):
    if alignment.aligned:
        coverage[alignment.iv] += 1
    i += 1

pyplot.plot(list(coverage[HTSeq.GenomicInterval("I", 200000, 500000, "+")]))
# pyplot.show()
'''

'''
file = "Homo_sapiens.GRCh38.113.abinitio.gtf"
gtf_file = HTSeq.GFF_Reader(file)

# this prints the features in reference genome
for feature in itertools.islice(gcf_file, 5):
    print("Feature: ", feature)
    print("Feature name: ", feature.name)
    print("Feature interval: ", feature.iv)
    print("Feature type: ", feature.type)
    print("Feature score: ", feature.score)
    print("Feature attributes: ", sorted(feature.attr.items()))
    print("="*100)
'''

'''
file = "Homo_sapiens.GRCh38.113.abinitio.gtf"
gtf_file = HTSeq.GFF_Reader(file)

exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)

# this loads in "exons" the names of the exons present in reference genome
for feature in gtf_file:
   if feature.type == "exon":
      exons[feature.iv] += feature.name

interval = HTSeq.GenomicInterval("III", 23850, 23950, ".")

iset = None

# this stores in iset the intersection of all the exons present in each step
for i, step_set in exons[interval].steps():
   if iset is None:
      iset = step_set.copy()
   else:
      iset.intersection_update(step_set)

print(iset)
'''

exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)

file = "Homo_sapiens.GRCh38.113.abinitio.gtf"
gtf_file = HTSeq.GFF_Reader(file)

# this loads in "exons" the names of the exons present in reference genome
for feature in gtf_file:
   if feature.type == "exon":
      exons[feature.iv] += feature.name

bam_reader = HTSeq.BAM_Reader("SRR1039508_sorted.bam")

counts = {}

for feature in gtf_file:
   if feature.type == "exon":
      counts[feature.name] = 0

for alignment in itertools.islice(bam_reader, 10000):
   if alignment.aligned:
      iset = None

      for iv2, step_set in exons[alignment.iv].steps():
          if iset is None:
             iset = step_set.copy()

          else:
             iset.intersection_update(step_set)

      if len(iset) == 1:
         counts[list(iset)[0]] += 1

for name in sorted(counts.keys()):
   print(name, counts[name])
