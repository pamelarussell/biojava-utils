package align;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import util.CommandLineParser;

import org.biojava.nbio.core.sequence.io.*;

/**
 * Perform pairwise alignment of two DNA sequences and write alignment information to file(s)
 * @author Pamela Russell
 *
 */
public class PairwiseAlign {
	
	// Logger
	private static Logger logger = LoggerFactory.getLogger(PairwiseAlign.class);
	
	// Two DNA sequences
	private DNASequence dna1;
	private DNASequence dna2;
	
	// Alignment algorithm parameters
	private SubstitutionMatrix<NucleotideCompound> subsMatrix = SubstitutionMatrixHelper.getNuc4_4();
	private SimpleGapPenalty gapPenalty = new SimpleGapPenalty();
	
	/**
	 * Construct with default alignment parameters
	 * @param seq1 DNA sequence 1
	 * @param seq2 DNA sequence 2
	 */
	public PairwiseAlign(DNASequence seq1, DNASequence seq2) {
		this.dna1 = seq1;
		this.dna2 = seq2;
	}
	
	/**
	 * Modify the gap open penalty
	 * @param g Positive gap open penalty
	 */
	public void setGapOpen(short g) {
		this.gapPenalty.setOpenPenalty(g);
	}
	
	/**
	 * Modify the gap extend penalty
	 * @param e Positive gap extend penalty
	 */
	public void setGapExtend(short e) {
		this.gapPenalty.setExtensionPenalty(e);
	}
	
	/**
	 * Returns the alignment of the two sequences
	 * @return Alignment
	 */
	public SequencePair<DNASequence, NucleotideCompound> getAlignment() {
		return Alignments.getPairwiseAlignment(dna1, dna2, PairwiseSequenceAlignerType.GLOBAL, this.gapPenalty, this.subsMatrix);
	}
	
	/**
	 * Main
	 * @param args
	 */
	public static void main(String[] args) {
				
		// Parse the command line
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-f1", "Fasta file 1 (contains one sequence only)", true);
		p.addStringArg("-f2", "Fasta file 2 (contains one sequence only)", true);
		p.parse(args);
		File fasta1 = new File(p.getStringArg("-f1"));
		File fasta2 = new File(p.getStringArg("-f2"));
		
		// Read the two sequences from the fasta files
		DNASequence seq1 = readSeqFromSingleSeqFasta(fasta1);
		DNASequence seq2 = readSeqFromSingleSeqFasta(fasta2);
		
		// Print the sequence lengths
		logger.info("Seq 1 length:\t" + seq1.getLength());
		logger.info("Seq 2 length:\t" + seq2.getLength());
		
		// Construct PairwiseAlign object
		PairwiseAlign pa = new PairwiseAlign(seq1, seq2);
		
		// Do the alignment
		SequencePair<DNASequence, NucleotideCompound> align = pa.getAlignment();
		
		// Analyze the alignment results
		List<AlignedSequence<DNASequence, NucleotideCompound>> alignments = align.getAlignedSequences();
		for(AlignedSequence<DNASequence, NucleotideCompound> a : align) {
			logger.info(Arrays.toString(a.getAlignmentFromSequence()));
			logger.info(Arrays.toString(a.getSequenceFromAlignment()));
		}
		
		logger.info("All done.");
		
	}
	
	/**
	 * Read a single DNA sequence from a fasta file
	 * Throws IllegalArgumentException if the file contains more than one sequence
	 * @param fasta Fasta file
	 * @return The single DNA sequence
	 */
	private static DNASequence readSeqFromSingleSeqFasta(File fasta) {
		Iterator<DNASequence> seqs;
		try {
			seqs = FastaReaderHelper.readFastaDNASequence(fasta).values().iterator();
			DNASequence rtrn = seqs.next();
			if(seqs.hasNext()) {
				throw new IllegalArgumentException("File " + fasta + " contains more than one sequence");
			}
			return rtrn;
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
			return null;
		}
	}
	
}
