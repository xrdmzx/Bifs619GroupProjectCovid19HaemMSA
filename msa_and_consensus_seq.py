from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo

# ORF1ab
# MSA using MUSCLE
in_file = "ORF1ab_regions_unique.fasta"
out_file = "ORF1ab_regions_aligned.fasta"
muscle_exe = "/Users/ronnyz/Documents/Fall2020/BIFS619/Group Project/heme_specific_datasets/muscle"
cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
stdout, stderr = cline(in_file)
align = AlignIO.read(out_file, "fasta")

# Consensus Sequence using 50% of sequences to call a position in the consensus sequence
alignment = AlignIO.read("ORF1ab_regions_aligned.fasta", 'fasta' )
summary_align = AlignInfo.SummaryInfo(alignment)
consensus_seq = summary_align.dumb_consensus(float(0.5))

# write consensus seqs to a text file
f = open( 'ORF1ab_consensus.seq.txt', 'w')
f.write(str(consensus_seq))
f.close()

##################################################################

# ORF3a
# MSA using MUSCLE
in_file = "ORF3a_regions_unique.fasta"
out_file = "ORF3a_regions_aligned.fasta"
muscle_exe = "/Users/ronnyz/Documents/Fall2020/BIFS619/Group Project/heme_specific_datasets/muscle"
cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
stdout, stderr = cline(in_file)
align = AlignIO.read(out_file, "fasta")

# Consensus Sequence using 50% of sequences to call a position in the consensus sequence
alignment = AlignIO.read("ORF3a_regions_aligned.fasta", 'fasta' )
summary_align = AlignInfo.SummaryInfo(alignment)
consensus_seq = summary_align.dumb_consensus(float(0.5))

# write consensus seqs to a text file
f = open( 'ORF3a_consensus.seq.txt', 'w')
f.write(str(consensus_seq))
f.close()

##################################################################

# ORF 8
# MSA using MUSCLE
in_file = "ORF8_regions_unique.fasta"
out_file = "ORF8_regions_aligned.fasta"
muscle_exe = "/Users/ronnyz/Documents/Fall2020/BIFS619/Group Project/heme_specific_datasets/muscle"
cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
stdout, stderr = cline(in_file)
align = AlignIO.read(out_file, "fasta")

# Consensus Sequence using 50% of sequences to call a position in the consensus sequence
alignment = AlignIO.read("ORF8_regions_aligned.fasta", 'fasta' )
summary_align = AlignInfo.SummaryInfo(alignment)
consensus_seq = summary_align.dumb_consensus(float(0.5))

# write consensus seqs to a text file
f = open( 'ORF8_consensus.seq.txt', 'w')
f.write(str(consensus_seq))
f.close()

##################################################################

# ORF10
# MSA using MUSCLE
in_file = "ORF10_regions_unique.fasta"
out_file = "ORF10_regions_aligned.fasta"
muscle_exe = "/Users/ronnyz/Documents/Fall2020/BIFS619/Group Project/heme_specific_datasets/muscle"
cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
stdout, stderr = cline(in_file)
align = AlignIO.read(out_file, "fasta")

# Consensus Sequence using 50% of sequences to call a position in the consensus sequence
alignment = AlignIO.read("ORF10_regions_aligned.fasta", 'fasta' )
summary_align = AlignInfo.SummaryInfo(alignment)
consensus_seq = summary_align.dumb_consensus(float(0.5))

# write consensus seqs to a text file
f = open( 'ORF10_consensus.seq.txt', 'w')
f.write(str(consensus_seq))
f.close()

##################################################################

# surface_glycoprotein
# MSA using MUSCLE
in_file = "surface_glycoprotein_regions_unique.fasta"
out_file = "surface_glycoprotein_regions_aligned.fasta"
muscle_exe = "/Users/ronnyz/Documents/Fall2020/BIFS619/Group Project/heme_specific_datasets/muscle"
cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
stdout, stderr = cline(in_file)
align = AlignIO.read(out_file, "fasta")

# Consensus Sequence using 50% of sequences to call a position in the consensus sequence
alignment = AlignIO.read("surface_glycoprotein_regions_aligned.fasta", 'fasta' )
summary_align = AlignInfo.SummaryInfo(alignment)
consensus_seq = summary_align.dumb_consensus(float(0.5))

# write consensus seqs to a text file
f = open( 'surface_glycoprotein_consensus.seq.txt', 'w')
f.write(str(consensus_seq))
f.close()
