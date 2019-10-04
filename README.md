# peakZCL

ChIP-Seq experiments are the most appropriate experimental choice for the study of gene
expression regulation. The critical step in the bioinformatic analysis of the obtained signals is the detection
of enriched regions or peaks. Although a great variety of algorithms are available, some challenges remain
to be addressed: a better definition of the peak morphology, a better quantification and statistical analysis,
and the development of tools for the comparison of the obtained results with different methods.

We propose a novel algorithm, peakZCL, based on the continuous wavelet transformof the ChIPSeq
signal. In order to provide a comprehensive performance evaluation we used two publicly available
experiments, EGR1 transcription factor and H3K4me3 histone modification ChIP-Seqs. In the case of
H3K4me3 signal, we generated a reference based on the manual annotation of the entire chromosome 16
provided by three experts. This resource, that allowed the comparison of peakZCL with other well-known
methods, can be of great interest for the research community. We achieved and excellent performance in
both case studies, outperforming most of the evaluated peak calling algorithms.
