
A modified GENOME coalescent simulator that outputs IBD segments. The original program can be found at http://www.sph.umich.edu/csg/liang/genome/
The simulator was modified to output IBD segments, defined as genomic regions for which pairs of samples share the same common ancestor. To use this feature, add the -ibd flag, followed by a number indicating the minimum length of the IBD segment to be output (e.g. "-ibd 1.0" for segments of 1cM or longer).
The function to output newick trees was also improved so that it is faster in this version.
IBD segments will be output using the format:
IBD   ID1   ID2   physicalFrom   physicalTo   geneticLength(cM)
Example:
./genome -pop 1 100 -N 1000 -pieces 10000 -ibd 1.0
The IBD segments are computed from the newick trees representing the ARG. An IBD segment is defined as a segment for which the most recent common ancestor for two haploid individuals remains unchanged. Tools to extract IBD using Newick trees in other simulators are available upon request.

Citations

- Original simulator: Liming Liang; Sebastian Zöllner; Goncalo R. Abecasis. "GENOME: a rapid coalescent-based whole genome simulator". Bioinformatics, 2007
- The IBD-enhanced simulator was developed in: P. F. Palamara, T. Lencz, A. Darvasi, I. Pe'er. "Length distributions of identity by descent reveal fine-scale demographic history". The American Journal of Human Genetics, 2012. 


Contacts

email Pier Palamara, ppalama AT hsph DOT harvard ONEMOREDOT edu for questions related to the -ibd flag.
Please contact the authors of the origina GENOME package for anything that is not related to the -ibd flag.

