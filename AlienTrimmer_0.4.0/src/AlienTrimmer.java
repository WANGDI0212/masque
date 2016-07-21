/*
  ###################################################################################
  # AlienTrimmer:  a  tool  to  quickly  and  accurately  trim  off  multiple short #
  # contaminant sequences from high-throughput sequencing reads                     #
  # [see version below]                                                             #
  # Copyright (C) 2014  Alexis Criscuolo                                            #
  #                                                                                 #
  # AlienTrimmer is free software;  you can redistribute it  and/or modify it under #
  # the terms of the  GNU General Public License  as published by the Free Software #
  # Foundation;  either version 2  of the  License,  or (at your option)  any later #
  # version.                                                                        #
  #                                                                                 #
  # AlienTrimmer is distributed in the hope that it will be useful, but WITHOUT ANY #
  # WARRANTY;  without even the implied warranty  of MERCHANTABILITY or FITNESS FOR #
  # A PARTICULAR PURPOSE. See the GNU General Public License for more details.      #
  #                                                                                 #
  # You should have received  a copy of the  GNU General Public License  along with #
  # this program; if not, write to the                                              #
  # Free Software Foundation, Inc.,                                                 #
  # 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA                         #
  #                                                                                 #
  # Contacts:                                                                       #
  # PF8 - Genotypage des Pathogenes et Sante Publique                               #
  # INSTITUT PASTEUR                                                                #
  # 28 rue du Dr Roux - 75724 Paris  (FRANCE)                                       #
  #                                                                                 #
  # GEM - Genomique Evolutive Microbienne                                           #
  # INSTITUT PASTEUR                                                                #
  # 28 rue du Dr Roux - 75724 Paris  (FRANCE)                                       #
  #                                                                                 #
  # alexis.criscuolo@pasteur.fr                                                     #
  # sylvain.brisse@pasteur.fr                                                       #
  ###################################################################################
*/
/*
  ####################################################################################################
  ### TECHNICAL NOTES ################################################################################
  ####################################################################################################
  # Different parameters:
    - s0:   first match within the read
    - i:    current index
    - cpt:  number of score entries > 0 up to index i
    - j = i+k-1

        s0          i         j
         |          |         |
  CCATCTGATCCCTGCGAGTCTCCGACTCAGGAGTTGCAGTGGTT  < read
         111111111 1234555554321                < 'score' array
  <--A--><----B----->
         <----------B'-------->

    =>  A = start0     B = i-start0+1     B' = j-start0+1

  # Computing 'start' such that start+1 is the 5' trimming index (initially, start := -1)
    If no k-mer match at index 'i', one have 4 criteria to set start := i (if start < i)
      (1) score[i] > 0
      (2) B >= A                  => i+1 >= 2*s0
      (3) cpt >= (A+B)/2          => 2*cpt >= i+1
      (4) i-start <= mismatch+1
    If a k-mer match occurs at index i, (1) is trivial and (4) remains;
    only (2) and (3) are modified to set start := j 
      (2') B' >= A                => j+1 >= 2*s0
      (3') cpt+k >= (A+B')/2      => 2*(cpt+k) >= j+1
  # Similar approach is used to compute 'end', the 3' trimming index (initially, end := read length)
  # Therefore, trimmed read is read.substring(start+1, end)
  ####################################################################################################
  ####################################################################################################
 */


import java.io.*;
import java.util.*;

public class AlienTrimmer {

    // constants
    static final String VERSION="0.4.0"; 
    static final byte B_1 = (byte) -1;
    static final byte B0 = (byte) 0;
    static final byte B1 = (byte) 1;
    static final byte B2 = (byte) 2;
    static final byte B3 = (byte) 3;
    static final byte B4 = (byte) 4;
    static final byte B5 = (byte) 5;
    static final byte B6 = (byte) 6;
    static final byte B7 = (byte) 7;
    static final byte B8 = (byte) 8;
    static final byte B9 = (byte) 9;
    static final byte B10 = (byte) 10;
    static final byte B15 = (byte) 15;
    static final char QINF = '!';
    static final long L0 = 0;
    
    // pre-specified alien sequence
    static final String   ALIEN0NAME = "Homopolymers";
    static final String[] ALIEN0 = { "AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGGGGGG", "TTTTTTTTTTTTTTTTTTTT" };
    static final String   ALIEN1NAME = "Dimers";
    static final String[] ALIEN1 = { "ACACACACACACACAC", "AGAGAGAGAGAGAGAG", "ATATATATATATATAT", 
				     "CGCGCGCGCGCGCGCG", "CTCTCTCTCTCTCTCT", "GTGTGTGTGTGTGTGT"  };
    static final String   ALIEN2NAME = "Trimers";
    static final String[] ALIEN2 = { "AACAACAACAACAACAAC", "AAGAAGAAGAAGAAGAAG", "AATAATAATAATAATAAT", 
				     "ACCACCACCACCACCACC", "ACGACGACGACGACGACG", "ACTACTACTACTACTACT",
				     "AGCAGCAGCAGCAGCAGC", "AGGAGGAGGAGGAGGAGG", "AGTAGTAGTAGTAGTAGT",
				     "ATCATCATCATCATCATC", "ATGATGATGATGATGATG", "ATTATTATTATTATTATT",
				     "CCGCCGCCGCCGCCGCCG", "CCTCCTCCTCCTCCTCCT", "CGGCGGCGGCGGCGGCGG",
				     "CGTCGTCGTCGTCGTCGT", "CTGCTGCTGCTGCTGCTG", "CTTCTTCTTCTTCTTCTT",
				     "GGTGGTGGTGGTGGTGGT", "GTTGTTGTTGTTGTTGTT" };
    static final String   ALIEN3NAME = "";
    static final String[] ALIEN3 = { };
    static final String   ALIEN4NAME = "";
    static final String[] ALIEN4 = { };
    static final String   ALIEN5NAME = "";
    static final String[] ALIEN5 = { };
    static final String   ALIEN6NAME = "";
    static final String[] ALIEN6 = { };
    static final String   ALIEN7NAME = "";
    static final String[] ALIEN7 = { };
    static final String   ALIEN8NAME = "";
    static final String[] ALIEN8 = { };
    static final String   ALIEN9NAME = "";
    static final String[] ALIEN9 = { };

    static void displayOptions() {
	System.out.println("");
	System.out.println(" AlienTrimmer v." + VERSION);
	System.out.println("");
	System.out.println(" USAGE:  AlienTrimmer  [options]");
	System.out.println("");
	System.out.println(" Fast trimming to filter out non-confident  nucleotides and alien oligo-");
	System.out.println(" nucleotide sequences (adaptors, primers) in both 5' and 3' read ends");
	System.out.println("");
	System.out.println(" OPTIONS:");
	System.out.println("  -i  <infile>  [single-ends] FASTQ formatted input file name");
	System.out.println("  -if <infile>  [paired-ends] FASTQ formatted input file name containing");
	System.out.println("                              forward (fwd) reads");
	System.out.println("  -ir <infile>  [paired-ends] FASTQ formatted input file name containing");
	System.out.println("                              reverse (rev) reads");
	System.out.println("  -o  <outfile> [single-ends] output file name");
	System.out.println("  -of <outfile> [paired-ends] output file name for trimmed fwd reads");
	System.out.println("  -or <outfile> [paired-ends] output file name for trimmed rev reads");
	System.out.println("  -os <outfile> [paired-ends] output  file  name for  remaining  trimmed");
	System.out.println("                              single (sgl) reads");
	System.out.println("  -c  [0-9]*    [single-ends] alien sequence  id(s) (see option -d)");
	System.out.println("  -d  [0-9]     displays alien sequences for the specified id:");
	if ( ALIEN0.length > 0 ) System.out.println("                 0: " + ALIEN0NAME);
	if ( ALIEN1.length > 0 ) System.out.println("                 1: " + ALIEN1NAME);
	if ( ALIEN2.length > 0 ) System.out.println("                 2: " + ALIEN2NAME);
	if ( ALIEN3.length > 0 ) System.out.println("                 3: " + ALIEN3NAME);
	if ( ALIEN4.length > 0 ) System.out.println("                 4: " + ALIEN4NAME);
	if ( ALIEN5.length > 0 ) System.out.println("                 5: " + ALIEN5NAME);
	if ( ALIEN6.length > 0 ) System.out.println("                 6: " + ALIEN6NAME);
	if ( ALIEN7.length > 0 ) System.out.println("                 7: " + ALIEN7NAME);
	if ( ALIEN8.length > 0 ) System.out.println("                 8: " + ALIEN8NAME);
	if ( ALIEN9.length > 0 ) System.out.println("                 9: " + ALIEN9NAME);
	System.out.println("  -c  <infile>  [single-ends] input file  name  containing  user-defined");
	System.out.println("                              alien  sequence(s) (one line per sequence)");
	System.out.println("  -cf [0-9]*    [paired-ends] same as -c for only fwd reads");
	System.out.println("  -cf <infile>  [paired-ends] same as -c for only fwd reads");
	System.out.println("  -cr [0-9]*    [paired-ends] same as -c for only rev reads");
	System.out.println("  -cr <infile>  [paired-ends] same as -c for only rev reads");
	System.out.println("  -k [5-15]     k value for k-mer decomposition;  must lie between 5 and");
	System.out.println("                15 (default: k=10)");
	System.out.println("  -m <int>      allowed mismatch value (default: m=k/2)");
	System.out.println("  -l <int>      minimum read  length to output;  all trimmed reads  with");
	System.out.println("                length below this value are filtered out (default: l=15)");
	System.out.println("  -q [0-40]     Phred quality score cutoff  to trim off low-quality read");
	System.out.println("                ends; must lie between 0 and 20 (default: q=20)");
	System.out.println("  -p [0-100]    minimum   allowed   percentage   of   correctly   called");
	System.out.println("                nucleotides  (i.e.  with Phred  quality score  character");
	System.out.println("                higher than q); all reads with a percentage of correctly");
	System.out.println("                called nucleotide lower than this value are filtered out");
	System.out.println("                (default: p=0)");
	System.out.println("  -v            displays trimming details during the whole process");
	System.out.println("");
	System.out.println(" EXAMPLES:");
	System.out.println(" [single-ends]");
	System.out.println("   AlienTrimmer -i reads.fq -o trim.fq -c 0 -l 30 -p 80");
	System.out.println("   AlienTrimmer -i reads.fq -c aliens.fa -k 9 -q 10");
	System.out.println(" [paired-ends]");
	System.out.println("   AlienTrimmer -if fwd.fq -ir rev.fq -c alien.fa -q 0 -p 75");
	System.out.println("   AlienTrimmer -if fwd.fq -ir rev.fq -cf alien.fwd.fa -cr alien.rev.fa");
	System.out.println("");
    }

    // options
    static File finfile;       // -i -if  fastq input file (fwd)
    static File rinfile;       // -ir     fastq input file (rev)
    static File foutfile;      // -o -of  fastq output file (fwd)
    static File routfile;      // -or     fastq output file (rev)
    static File soutfile;      // -os     fastq output file (sgl)
    static String falien;      // -c -cf  contaminant library name, or user file name (fwd)
    static String ralien;      // -cr     contaminant library name, or user file name (rev)
    static byte k;             // -k      value for contaminant k-mer decomposition
    static byte mismatch;      // -m      nb of allowed mismatches
    static char qmin;          // -q      minimum allowed Q score
    static short prop;         // -p      minimum percentage of correct nucleotides (i.e. >= qmin)
    static int minLgt;         // -l      min read length
    static byte verbose;       // -v      to output trimming details
    
    // io
    static BufferedReader fin, rin;
    static BufferedWriter fout, rout, sout;

    // data
    static int phred;
    static byte paired;
    static BitSet bsfkmer, bsrkmer;
    static int[] ifkmer, irkmer;
    static long[] lfkmer, lrkmer;
    static String fid1, fseq, fid2, fqsc, rid1, rseq, rid2, rqsc;
    static byte[] score;

    // stuffs
    static byte cpt, b, bq, k_1, k2;
    static short i, j, x, lgt, l, start0, start, end0, end, n, scount;
    static int qval, id, o, lfkmerLgt, lrkmerLgt;
    static int pcpt, ftcpt, rtcpt, frcpt, rrcpt;
    static long t, cur;
    static int imsk, iptn;
    static long lmsk, lptn, nptn, lkm, lpt;
    static String line, seq;
    static ArrayList<String> contaminant;
    static TreeSet<Long> lts;
    static StringBuffer sb;

    public static void main(String[] args) throws IOException {

	if ( args.length < 2 ) { displayOptions(); System.exit(0); }

	//#####################
	//#####################
	//## reading options ##
	//#####################
	//#####################
	paired = B0;
	finfile = new File("no.file"); rinfile = new File("no.file");
	foutfile = new File("no.file"); routfile = new File("no.file"); soutfile = new File("no.file"); 
	falien = "no.alien"; ralien = "no.alien";
	k = B10;
	mismatch = B_1;
	minLgt = 15;
	qval = 20; //qmin = QINF;
	prop = B0;
	verbose = B0;
	o = -1;
	while ( ++o < args.length ) {
	    //### -i  fastq input file (fwd) #################################################################################
	    if ( args[o].equals("-i") ) {
		finfile = new File(args[++o]);
		if ( ! finfile.exists() ) { System.out.println("  input file does not exists (option -i)"); System.exit(0); }
		continue;
	    }
	    //### -if  fastq input file (fwd) ################################################################################
	    if ( args[o].equals("-if") ) {
		finfile = new File(args[++o]);
		if ( ! finfile.exists() ) { System.out.println("  input file does not exists (option -if)"); System.exit(0); }
		continue;
	    }
	    //### -ir  fastq input file (rev) ################################################################################
	    if ( args[o].equals("-ir") ) {
		rinfile = new File(args[++o]);
		if ( ! rinfile.exists() ) { System.out.println("  input file does not exists (option -ir)"); System.exit(0); }
		paired = B1;
		continue;
	    }
	    //### -o  fastq output file (fwd) ################################################################################
	    if ( args[o].equals("-o") ) { foutfile = new File(args[++o]); continue; }
	    //### -of  fastq output file (fwd) ###############################################################################
	    if ( args[o].equals("-of") ) { foutfile = new File(args[++o]); continue; }
	    //### -or  fastq output file (rev) ###############################################################################
	    if ( args[o].equals("-or") ) { routfile = new File(args[++o]); paired = B1; continue; }
	    //### -os  fastq output file (sgl) ###############################################################################
	    if ( args[o].equals("-os") ) { soutfile = new File(args[++o]); paired = B1; continue; }
	    //### -c  contaminant id, or user file name (fwd) ################################################################
	    if ( args[o].equals("-c") ) { falien = args[++o]; continue; }
	    //### -cf  contaminant id, or user file name (fwd) ###############################################################
	    if ( args[o].equals("-cf") ) { falien = args[++o]; continue; }
	    //### -cr  contaminant id, or user file name (rev) ###############################################################
	    if ( args[o].equals("-cr") ) { ralien = args[++o]; paired = B1; continue; }
	    //### -k  value for k-mer decomposition ##########################################################################
	    if ( args[o].equals("-k") ) { 
		try { k = Byte.parseByte(args[++o]); } 
		catch ( NumberFormatException e ) { System.out.println("  incorrect K-mer value (option -k)"); System.exit(0); }
		if ( (k < B5) || (k > B15) ) { System.out.println("  K-mer should be set between 5 and 15 (option -k)"); System.exit(1); }
		continue;
	    }
	    //### -m  nb of allowed mismatches ###############################################################################
	    if ( args[o].equals("-m") ) { 
		try { mismatch = Byte.parseByte(args[++o]); } 
		catch ( NumberFormatException e ) { System.out.println("  incorrect mismatch value (option -m)"); System.exit(0); }
		if ( (mismatch < B0) || (mismatch > B15) ) { System.out.println("  mismatch value should be set between 0 and 15 (option -m)"); System.exit(1); }
		continue;
	    }
	    //### -l  min read length ########################################################################################
	    if ( args[o].equals("-l") ) { 
		try { minLgt = Integer.parseInt(args[++o]); } 
		catch ( NumberFormatException e ) { System.out.println("  incorrect min read length value (option -l)"); System.exit(0); }
		if ( minLgt < 0 ) { System.out.println("  min read length should be positive (option -l)"); System.exit(1); }
		continue;
	    }
	    //### -q  minimum allowed Q score #################################################################################
	    if ( args[o].equals("-q") ) { 
		try { qval = Integer.parseInt(args[++o]); } 
		catch ( NumberFormatException e ) { System.out.println("  incorrect Phred quality score cutoff (option -q)"); System.exit(0); } 
		if ( (qval < 0) || (qval > 41) ) { System.out.println("  Phred quality score cutoff should be set between 0 and 40 (option -q)"); System.exit(0); } 
		continue;
	    }
	    //### -p   minimum percentage of correct nucleotides (i.e. >= qmin) ###############################################
	    if ( args[o].equals("-p") ) { 
		try { prop = (short) Double.parseDouble(args[++o]); }
		catch ( NumberFormatException e ) { System.out.println("  incorrect percentage value (option -p)"); System.exit(0); }
		if ( prop < 0 ) prop = B0;
		if ( prop > 100 ) { System.out.println("  percentage should be set between 0 and 100 (option -p)"); System.exit(1); }
		continue;
	    }
	    //### -v  verbose ################################################################################################
	    if ( args[o].equals("-v") ) { verbose = B1; continue; }
	    //### -d  displaying contaminants ################################################################################
	    if ( args[o].equals("-d") ) {
		try { id = Integer.parseInt(args[++o]); } 
		catch ( NumberFormatException e ) { System.out.println("  incorrect alien sequence id (option -d)"); System.exit(1); }
		contaminant = new ArrayList<String>();
		switch ( id ) {
		case 0: Collections.addAll(contaminant, ALIEN0); break;
		case 1: Collections.addAll(contaminant, ALIEN1); break;
		case 2: Collections.addAll(contaminant, ALIEN2); break;
		case 3: Collections.addAll(contaminant, ALIEN3); break;
		case 4: Collections.addAll(contaminant, ALIEN4); break;
		case 5: Collections.addAll(contaminant, ALIEN5); break;
		case 6: Collections.addAll(contaminant, ALIEN6); break;
		case 7: Collections.addAll(contaminant, ALIEN7); break;
		case 8: Collections.addAll(contaminant, ALIEN8); break;
		case 9: Collections.addAll(contaminant, ALIEN9); break;
		}
		id = -1; while ( ++id < contaminant.size() ) System.out.println(contaminant.get(id));
		System.exit(0);
	    }
	}
	//### testing default variables ##############################################################################
	if ( finfile.toString().equals("no.file") ) { System.out.println("  no input file"); System.exit(1); }
	if ( foutfile.toString().equals("no.file") ) foutfile = new File(finfile + ".at.fq");
	if ( paired != B0 ) {
	    if ( rinfile.toString().equals("no.file") ) { System.out.println("  no input file (option -ir)"); System.exit(1); }
	    if ( routfile.toString().equals("no.file") ) routfile = new File(rinfile + ".at.fq");
	    if ( soutfile.toString().equals("no.file") ) soutfile = new File(finfile + ".at.sgl.fq");
	}
	if ( mismatch < B0 ) { mismatch = k; ++mismatch; mismatch /= B2; }
	//### detecting Phred encoding ##############################################################################
	phred = 0; fin = new BufferedReader(new FileReader(finfile)); 
	cpt = B0; 
	while ( phred == 0 ) {		
	    //## reading fastq ##
	    try { line = fin.readLine().trim(); } catch ( NullPointerException e ) { fin.close(); break; }
	    switch ( ++cpt ) { case B1: case B2: case B3: continue; case B4: fqsc = line; cpt = B0; break; }
	    i = (short) fqsc.length(); while ( --i >= 0 ) { if ( fqsc.charAt(i) <= '9' ) { phred = 33; break; } if ( fqsc.charAt(i) >= 'L' ) { phred = 64; break; } }
	}
	qmin = (char) (phred+qval);
	if ( qval <= 0 ) prop = B0;

	System.out.print("AlienTrimmer main options: -k " + k + " -l " + minLgt + " -m " + mismatch + " -q " + qval + " -p " + prop + " (Phred+" + phred + ")  /  "); 


	//####################
	//####################
	//## init variables ##
	//####################
	//####################
	t = System.currentTimeMillis();
	imsk = 0; ++imsk; imsk <<= B2*k; --imsk;
	lmsk = L0; ++lmsk; lmsk <<= B4*k; --lmsk; 
	k_1 = k; --k_1;
	k2 = k; k2 *= B2;
	++mismatch;
	

	//###################################################
	//###################################################
	//## single-ends reads                             ##
	//###################################################
	//###################################################
	if ( paired == B0 ) {
	    //##################################
	    //## computing contaminant k-mers ##
	    //##################################
	    contaminant = new ArrayList<String>();
	    if ( (new File(falien)).exists() ) {
		fin = new BufferedReader(new FileReader(new File(falien)));
		while ( true ) {
		    try { line = fin.readLine().trim(); } catch ( NullPointerException e ) { fin.close(); break; }
		    if ( (line.length() == 0) || line.startsWith("#") || line.startsWith("%") || line.startsWith(">")) continue;
		    contaminant.add( line );
		}
	    }
	    else {
		o = falien.length();
		while ( --o >= 0 ) {
		    try { id = Integer.parseInt("" + falien.charAt(o)); }
		    catch ( NumberFormatException e ) { System.out.println("  incorrect alien sequence id (option -c)"); System.exit(0); }
		    switch ( id ) {
		    case 0: Collections.addAll(contaminant, ALIEN0); break; case 1: Collections.addAll(contaminant, ALIEN1); break; 
		    case 2: Collections.addAll(contaminant, ALIEN2); break; case 3: Collections.addAll(contaminant, ALIEN3); break;
		    case 4: Collections.addAll(contaminant, ALIEN4); break; case 5: Collections.addAll(contaminant, ALIEN5); break;
		    case 6: Collections.addAll(contaminant, ALIEN6); break; case 7: Collections.addAll(contaminant, ALIEN7); break;
		    case 8: Collections.addAll(contaminant, ALIEN8); break; case 9: Collections.addAll(contaminant, ALIEN9); break;
		    }
		}
	    }
	    if ( contaminant.size() == 0 ) { System.out.println("  no alien sequence (option -c)"); System.exit(0); }
	    System.out.print(contaminant.size() + " alien sequence(s)  /  ");
	    //## verifying each contaminant nucleotides ##
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		sb = new StringBuffer(contaminant.get(o)); 
		if ( (i=(short)sb.length()) < k ) { contaminant.remove(o); continue; }
		while ( --i >= 0 ) {
		    switch ( sb.charAt(i) ) {
		    case 'A': case 'a': case 'C': case 'c': case 'G': case 'g': case 'T': case 't': continue;
		    case 'U': case 'u': sb.setCharAt(i, 'T'); contaminant.set(o, sb.toString()); break;
		    case 'M': case 'm': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'R': case 'r': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'W': case 'w': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'S': case 's': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'Y': case 'y': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'K': case 'k': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'B': case 'b': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'D': case 'd': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'H': case 'h': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'V': case 'v': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'N': case 'n': case 'X': case 'x': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); 
			sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    default: System.out.println("the following alien sequence contains at least one incorrect nucleotide (" + sb.charAt(i) + ") : " + sb.toString()); System.exit(1);
		    }
		}
	    }
	    //## reverse-complementing ##
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		sb = new StringBuffer((lgt=(short)(seq=contaminant.get(o)).length()));
		while ( --lgt >= 0 ) 
		    switch ( seq.charAt(lgt) ) {
		    case 'a': case 'A': sb = sb.append('T'); break; case 'c': case 'C': sb = sb.append('G'); break; 
		    case 'g': case 'G': sb = sb.append('C'); break; case 't': case 'T': sb = sb.append('A'); break; 
		    }
		contaminant.add( sb.toString() );
	    }
	    //## k-mering ##
	    bsfkmer = new BitSet(imsk+1); lts = new TreeSet<Long>();
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		lgt = (short) (seq=contaminant.get(o).toUpperCase()).length(); 
		iptn = B0; lptn = L0; j = B_1;
		while ( ++j < k_1 ) { 
		    iptn <<= B2; lptn <<= B4; 
		    switch ( seq.charAt(j) ) { 
		    case 'A':             lptn |= B1; break; case 'C': iptn |= B1; lptn |= B2; break; 
		    case 'G': iptn |= B2; lptn |= B4; break; case 'T': iptn |= B3; lptn |= B8; break; 
		    }		    
		}
		--j;
		while ( ++j < lgt ) {
		    iptn <<= B2; lptn <<= B4; 
		    switch ( seq.charAt(j) ) { 
		    case 'A':             lptn |= B1; break; case 'C': iptn |= B1; lptn |= B2; break; 
		    case 'G': iptn |= B2; lptn |= B4; break; case 'T': iptn |= B3; lptn |= B8; break; 
		    }
		    bsfkmer.set((iptn &= imsk)); lts.add( new Long((lptn &= lmsk)) );
		}
	    }
	    lfkmer = new long[(lfkmerLgt=(o=lts.size()))]; while ( --o >= 0 ) lfkmer[o] = lts.pollLast().longValue(); 
	    System.out.println(lfkmer.length + " k-mers (k=" + k + ")"); 
	    lts = null; contaminant = null; sb = null;

	    //####################
	    //## trimming reads ##
	    //####################
	    fin = new BufferedReader(new FileReader(finfile)); fout = new BufferedWriter(new FileWriter(foutfile));
	    pcpt = 0; ftcpt = 0; frcpt = 0;
	    cpt = B0; score = new byte[51];
	    while ( true ) {		
		//## reading fastq ##
		try { line = fin.readLine().trim(); } catch ( NullPointerException e ) { fin.close(); break; }
		switch ( ++cpt ) { case B1: fid1 = line; continue; case B2: fseq = line; continue; case B3: fid2 = line; continue; case B4: fqsc = line; cpt = B0; break; }
		
		//## displaying info ##
		if ( ++pcpt % 1000000 == 0 ) 
		    System.out.println("[" + String.format(Locale.US, "%02d", new Long((cur=(System.currentTimeMillis()-t)/1000)/60))
				       + ":" + (line="0"+(cur%60)).substring(line.length()-2) + "]" 
				       + String.format(Locale.US, "%,12d", new Integer(pcpt)) + " reads processed:" 
				       + String.format(Locale.US, "%,12d", new Integer(ftcpt)) + " trimmed" 
				       + String.format(Locale.US, "%,12d", new Integer(frcpt)) + " removed");
		
		//## matching k-mers... ##
		if ( (lgt=(short)fseq.length()) < minLgt ) { ++frcpt; continue; } // too small read
		if ( lgt >= score.length ) { score = new byte[++lgt]; --lgt; } else Arrays.fill(score, B0);
		end = lgt; start = B_1; start0 = B_1; scount = 0; iptn = B0; i = B_1; j = B_1; b = B0; bq = B0; 
		if ( (fseq.indexOf('N') != -1) /*|| (fseq.indexOf('n') != -1)*/ ) {
		    //## ... with degenerate nucleotides ##
		    lptn = L0; nptn = L0; n = B0;
		    while ( ++j < k_1 ) { 
			iptn <<= B2; lptn <<= B4; nptn <<= B4; 
			switch ( fseq.charAt(j) ) { 
			case 'A': case 'a':             lptn |= B1; nptn |= B1;  break; 
			case 'C': case 'c': iptn |= B1; lptn |= B2; nptn |= B2;  break; 
			case 'G': case 'g': iptn |= B2; lptn |= B4; nptn |= B4;  break; 
			case 'T': case 't': iptn |= B3; lptn |= B8; nptn |= B8;  break; 
			default:  n = k;                            nptn |= B15; 
			} 
			--n;
		    }
		    while ( j < lgt ) {
			b = (score[++i] += b); bq = ( fqsc.charAt(i) < qmin ) ? B1 : B0;
			iptn <<= B2; lptn <<= B4; nptn <<= B4; 
			switch ( fseq.charAt(j) ) { 
			case 'A': case 'a':             lptn |= B1; nptn |= B1;  break; 
			case 'C': case 'c': iptn |= B1; lptn |= B2; nptn |= B2;  break; 
			case 'G': case 'g': iptn |= B2; lptn |= B4; nptn |= B4;  break; 
			case 'T': case 't': iptn |= B3; lptn |= B8; nptn |= B8;  break; 
			default:  n = k;                            nptn |= B15; 
			} 
			switch ( (((--n) < B0) ? B0 : B1) ) {
			case B0:
			    if ( bsfkmer.get((iptn &= imsk)) ) { 
				b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
				scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
				//start = ( /*(i-start <= k2) &&*/ (j >= start0) && (scount+k2 >= j) ) ? j : ++start; --start;
				start = ( (j >= start0) && (scount+k2 >= j) ) 
				    ? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				    : ++start; --start;
				continue; 
			    }
			    break;
			default:
			    id = Arrays.binarySearch(lfkmer, (lptn &= lmsk)); id = -(++id); --id;
			    o = Arrays.binarySearch(lfkmer, (nptn &= lmsk)); o = -(++o); 
			    while ( (++id < o) && ((lkm=lfkmer[id]) != (lkm |= (lpt=lptn))) ) {}
			    if ( id != o ) { 
				b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
				scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
				//start = ( /*(i-start <= k2) &&*/ (j >= start0) && (scount+k2 >= j) ) ? j : ++start; --start;
				start = ( (j >= start0) && (scount+k2 >= j) ) 
				    ? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				    : ++start; --start;
				continue; 
			    }
			    break;
			}
			b = (score[i]=((b < B0) ? B0 : b)); ++j; 
			if ( (score[i] += bq) == B0 ) continue;
			scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
			start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		    }
		}
		else { 
		    //## ... with no degenerate nucleotide ##
		    while ( ++j < k_1 ) { 
			iptn <<= B2; 
			switch ( fseq.charAt(j) ) { case 'C': case 'c': iptn |= B1; break; case 'G': case 'g': iptn |= B2; break; case 'T': case 't': iptn |= B3; break; } 
		    }
		    while ( j < lgt ) {
			//System.out.println(i + " " + j + " " + (start0/3) + " " + start + " " + (scount/2) + " " + end);
			b = (score[++i] += b); bq = ( fqsc.charAt(i) < qmin ) ? B1 : B0;
			iptn <<= B2; switch ( fseq.charAt(j) ) { case 'C': case 'c': iptn |= B1; break; case 'G': case 'g': iptn |= B2; break; case 'T': case 't': iptn |= B3; break; } 
			if ( bsfkmer.get((iptn &= imsk)) ) { 
			    b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
			    scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
			    //start = ( /*(i-start <= k2) &&*/ (j >= start0) && (scount+k2 >= j) ) ? j : ++start; --start;
			    start = ( (j >= start0) && (scount+k2 >= j) ) 
				? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				: ++start; --start;
			    continue; 
			}
			b = (score[i]=((b < B0) ? B0 : b)); ++j; 
			if ( (score[i] += bq) == B0 ) continue;
			scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
			start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		    }
		}
		while ( ++i < lgt ) { 
		    b = (score[i] += b); bq = ( fqsc.charAt(i) < qmin ) ? B1 : B0; b = (score[i]=((b < B0) ? B0 : b)); 
		    if ( (score[i] += bq) == B0 ) continue;
		    scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
		    start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		}
		start = ( (start0 > B0) && (start0 >= lgt-start) ) ? B_1 : start;

		// here :
		//   start0/3:    index of the first non-zero score value
		//   start+1:   5' trimming index
		//   scount/2:  number of non-zero score values up to start index
		//   end:       index of the last non-zero score value
		// therefore, one only have to look for 3' trimming index
		if ( (start0 >= B0) && (lgt-start-B1 >= minLgt) ) {
		    end0 = end; end = lgt; 
		    j = end0; i = j; i -= k; ++i; 
		    ++end0; end0 *= B2;
		    scount = B2;
		    end = ( (end > j) && (end-j <= mismatch) && (end0-j >= lgt) && (scount >= lgt-j) ) ? j : end;
		    start0 /= B2; o = ( start0 > start ) ? start0 : start; 
		    //System.out.println("" + start0 + " " + end + " " + lgt);
		    while ( (x=(--j)) >= o ) {
			--i; 
			if ( (b=score[j]) == B0 ) continue;
			scount += B2;
			/*
			if ( (b >= B2) && ((--b) == score[++x]) && (score[i] != B0) ) 
			    end = ( (end0-i >= lgt) && (scount+k2 >= lgt-i) ) 
			    ? (( end == lgt ) ? i : (( end-j <= mismatch ) ? i : end)) 
				: end; 
			*/
			end = ( (i >= B0) && (b >= B2) && ((--b) == score[++x]) && (score[i] != B0) ) 
			      ? ( (end0-i >= lgt) && (scount+k2 >= lgt-i) ) 
			        ? (( end == lgt ) ? i : (( end-j <= mismatch ) ? i : end)) 
			        : end
			      : end ;
			end = ( (end > j) && (end-j <= mismatch) && (end0-j >= lgt) && (scount >= lgt-j) ) ? j : end;
			//System.out.println(i + " " + j + " " + x + "  " + (scount/2) + "  " + (end0/2-1) + " " + end);
		    }
		    end = ( start < end ) ? --end : start; ++end;
		}
		else end = lgt;

		// here the trimmed read is read.substring(++start, end)
		l = end; l -= start; --l;

		// verifying composition of remaining incorrect nucleotides
		scount = l; if ( prop != B0 ) { i = start; while ( ++i < end ) scount -= ( score[i] == B1 ) ? B1 : B0; }

		//## writing trimmed read ##
		if ( (l < minLgt) || (100*scount < prop*l) ) ++frcpt;
		else {
		    if ( (start == B_1) && (end == lgt) ) {
			fout.write(fid1); fout.newLine(); fout.write(fseq); fout.newLine();
			fout.write(fid2); fout.newLine(); fout.write(fqsc); fout.newLine();
		    }
		    else {
			++ftcpt;
			fout.write(fid1); fout.newLine(); fout.write(fseq.substring(++start,end)); fout.newLine();
			fout.write(fid2); fout.newLine(); fout.write(fqsc.substring(start--,end)); fout.newLine();
		    }
		}

		//## displaying results ##
		if ( (verbose == B0) || ((start == B_1) && (end == lgt)) ) continue;
		System.out.println(fid1); System.out.println(fseq); //if ( qmin != QINF ) { System.out.println("+"); System.out.println(fqsc); }
		o = B_1; while ( ++o < lgt ) System.out.print( (( (b=score[o]) == B0 ) ? " " : ( b < B10 ) ? b : ( b == B10 ) ? "0" : (""+(char)(54+b))) ); System.out.println("");
		o = B_1; while ( ++o <= start ) System.out.print(">"); --o; while ( ++o < end ) System.out.print(" ");  --o; while ( ++o < lgt ) System.out.print("<"); 
		System.out.println(""); System.out.println("");
	    }	
	    fout.close();
	    System.out.println("[" + String.format(Locale.US, "%02d", new Long((cur=(System.currentTimeMillis()-t)/1000)/60))
			       + ":" + (line="0"+(cur%60)).substring(line.length()-2) + "]" 
			       + String.format(Locale.US, "%,12d", new Integer(pcpt)) + " reads processed:" 
			       + String.format(Locale.US, "%,12d", new Integer(ftcpt)) + " trimmed" 
			       + String.format(Locale.US, "%,12d", new Integer(frcpt)) + " removed");
	}





	//###################################################
	//###################################################
	//## paired-ends reads                             ##
	//###################################################
	//###################################################
	if ( paired != B0 ) {
	    //########################################
	    //## computing contaminant K-mers (fwd) ##
	    //########################################
	    contaminant = new ArrayList<String>();
	    if ( (new File(falien)).exists() ) {
		fin = new BufferedReader(new FileReader(new File(falien)));
		while ( true ) {
		    try { line = fin.readLine().trim(); } catch ( NullPointerException e ) { fin.close(); break; }
		    if ( (line.length() == 0) || line.startsWith("#") || line.startsWith("%")  || line.startsWith(">")) continue;
		    contaminant.add( line );
		}
	    }
	    else {
		o = falien.length();
		while ( --o >= 0 ) {
		    try { id = Integer.parseInt("" + falien.charAt(o)); }
		    catch ( NumberFormatException e ) { System.out.println("  incorrect alien sequence id (option -cf)"); System.exit(0); }
		    switch ( id ) {
		    case 0: Collections.addAll(contaminant, ALIEN0); break;
		    case 1: Collections.addAll(contaminant, ALIEN1); break;
		    case 2: Collections.addAll(contaminant, ALIEN2); break;
		    case 3: Collections.addAll(contaminant, ALIEN3); break;
		    case 4: Collections.addAll(contaminant, ALIEN4); break;
		    case 5: Collections.addAll(contaminant, ALIEN5); break;
		    case 6: Collections.addAll(contaminant, ALIEN6); break;
		    case 7: Collections.addAll(contaminant, ALIEN7); break;
		    case 8: Collections.addAll(contaminant, ALIEN8); break;
		    case 9: Collections.addAll(contaminant, ALIEN9); break;
		    }
		}
	    }
	    if ( contaminant.size() == 0 ) { System.out.println("  no alien sequence  (option -cf)"); System.exit(0); }
	    System.out.print(contaminant.size() + " fwd alien sequence(s)  /  ");
	    //## verifying each contaminant nucleotides ##
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		sb = new StringBuffer(contaminant.get(o)); 
		if ( (i=(short)sb.length()) < k ) { contaminant.remove(o); continue; }
		while ( --i >= 0 ) {
		    switch ( sb.charAt(i) ) {
		    case 'A': case 'a': case 'C': case 'c': case 'G': case 'g': case 'T': case 't': continue;
		    case 'U': case 'u': sb.setCharAt(i, 'T'); contaminant.set(o, sb.toString()); break;
		    case 'M': case 'm': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'R': case 'r': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'W': case 'w': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'S': case 's': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'Y': case 'y': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'K': case 'k': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'B': case 'b': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'D': case 'd': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'H': case 'h': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'V': case 'v': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'N': case 'n': case 'X': case 'x': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); 
			sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    default: System.out.println("the following alien sequence contains at least one incorrect nucleotide (" + sb.charAt(i) + ") : " + sb.toString()); System.exit(1);
		    }
		}
	    }
	    //## reverse-complementing ##
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		sb = new StringBuffer((lgt=(short)(seq=contaminant.get(o)).length()));
		while ( --lgt >= 0 ) 
		    switch ( seq.charAt(lgt) ) {
		    case 'a': case 'A': sb = sb.append('T'); break; case 'c': case 'C': sb = sb.append('G'); break; 
		    case 'g': case 'G': sb = sb.append('C'); break; case 't': case 'T': sb = sb.append('A'); break; 
		    }
		contaminant.add( sb.toString() );
	    }
	    //## k-mering ##
	    bsfkmer = new BitSet(imsk+1); lts = new TreeSet<Long>();
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		lgt = (short) (seq=contaminant.get(o).toUpperCase()).length(); 
		iptn = B0; lptn = L0; j = B_1;
		while ( ++j < k_1 ) { 
		    iptn <<= B2; lptn <<= B4; 
		    switch ( seq.charAt(j) ) { 
		    case 'A':             lptn |= B1; break; case 'C': iptn |= B1; lptn |= B2; break; 
		    case 'G': iptn |= B2; lptn |= B4; break; case 'T': iptn |= B3; lptn |= B8; break; 
		    }		    
		}
		--j;
		while ( ++j < lgt ) {
		    iptn <<= B2; lptn <<= B4; 
		    switch ( seq.charAt(j) ) { 
		    case 'A':             lptn |= B1; break; case 'C': iptn |= B1; lptn |= B2; break; 
		    case 'G': iptn |= B2; lptn |= B4; break; case 'T': iptn |= B3; lptn |= B8; break; 
		    }
		    bsfkmer.set((iptn &= imsk)); lts.add( new Long((lptn &= lmsk)) );
		}
	    }
	    lfkmer = new long[(lfkmerLgt=(o=lts.size()))]; while ( --o >= 0 ) lfkmer[o] = lts.pollLast().longValue(); 
	    System.out.print(lfkmer.length + " fwd k-mers (k=" + k + ")  /  "); 

	    //########################################
	    //## computing contaminant K-mers (rev) ##
	    //########################################
	    contaminant = new ArrayList<String>();
	    if ( ralien.equals("no.alien") ) { try { id = Integer.parseInt(ralien); } catch ( NumberFormatException e ) { ralien = falien; } }
	    if ( (new File(ralien)).exists() ) {
		rin = new BufferedReader(new FileReader(new File(ralien)));
		while ( true ) {
		    try { line = rin.readLine().trim(); } catch ( NullPointerException e ) { rin.close(); break; }
		    if ( (line.length() == 0) || line.startsWith("#") || line.startsWith("%")  || line.startsWith(">")) continue;
		    contaminant.add( line );
		}
	    }
	    else {
		o = ralien.length();
		while ( --o >= 0 ) {
		    try { id = Integer.parseInt("" + ralien.charAt(o)); }
		    catch ( NumberFormatException e ) { System.out.println("  incorrect alien sequence id (option -cr)"); System.exit(0); }
		    switch ( id ) {
		    case 0: Collections.addAll(contaminant, ALIEN0); break;
		    case 1: Collections.addAll(contaminant, ALIEN1); break;
		    case 2: Collections.addAll(contaminant, ALIEN2); break;
		    case 3: Collections.addAll(contaminant, ALIEN3); break;
		    case 4: Collections.addAll(contaminant, ALIEN4); break;
		    case 5: Collections.addAll(contaminant, ALIEN5); break;
		    case 6: Collections.addAll(contaminant, ALIEN6); break;
		    case 7: Collections.addAll(contaminant, ALIEN7); break;
		    case 8: Collections.addAll(contaminant, ALIEN8); break;
		    case 9: Collections.addAll(contaminant, ALIEN9); break;
		    }
		}
	    }
	    if ( contaminant.size() == 0 ) { System.out.println("  no alien sequence (option -cr)"); System.exit(0); }
	    System.out.print(contaminant.size() + " rev alien sequence(s)  /  ");
	    //## verifying each contaminant nucleotides ##
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		sb = new StringBuffer(contaminant.get(o)); 
		if ( (i=(short)sb.length()) < k ) { contaminant.remove(o); continue; }
		while ( --i >= 0 ) {
		    switch ( sb.charAt(i) ) {
		    case 'A': case 'a': case 'C': case 'c': case 'G': case 'g': case 'T': case 't': continue;
		    case 'U': case 'u': sb.setCharAt(i, 'T'); contaminant.set(o, sb.toString()); break;
		    case 'M': case 'm': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'R': case 'r': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'W': case 'w': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'S': case 's': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'Y': case 'y': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'K': case 'k': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'B': case 'b': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'D': case 'd': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'H': case 'h': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'V': case 'v': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    case 'N': case 'n': case 'X': case 'x': 
			contaminant.remove(o); 
			sb.setCharAt(i, 'A'); contaminant.add(sb.toString()); sb.setCharAt(i, 'C'); contaminant.add(sb.toString()); 
			sb.setCharAt(i, 'G'); contaminant.add(sb.toString()); sb.setCharAt(i, 'T'); contaminant.add(sb.toString()); 
			o = contaminant.size(); i = 0; break;
		    default: System.out.println("the following alien sequence contains at least one incorrect nucleotide (" + sb.charAt(i) + ") : " + sb.toString()); System.exit(1);
		    }
		}
	    }
	    //## reverse-complementing ##
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		sb = new StringBuffer((lgt=(short)(seq=contaminant.get(o)).length()));
		while ( --lgt >= 0 ) 
		    switch ( seq.charAt(lgt) ) {
		    case 'a': case 'A': sb = sb.append('T'); break; case 'c': case 'C': sb = sb.append('G'); break; 
		    case 'g': case 'G': sb = sb.append('C'); break; case 't': case 'T': sb = sb.append('A'); break; 
		    }
		contaminant.add( sb.toString() );
	    }
	    //## k-mering ##
	    bsrkmer = new BitSet(imsk+1); lts = new TreeSet<Long>();
	    o = contaminant.size();
	    while ( --o >= 0 ) {
		lgt = (short) (seq=contaminant.get(o).toUpperCase()).length(); 
		iptn = B0; lptn = L0; j = B_1;
		while ( ++j < k_1 ) { 
		    iptn <<= B2; lptn <<= B4; 
		    switch ( seq.charAt(j) ) { 
		    case 'A':             lptn |= B1; break; case 'C': iptn |= B1; lptn |= B2; break; 
		    case 'G': iptn |= B2; lptn |= B4; break; case 'T': iptn |= B3; lptn |= B8; break; 
		    }		    
		}
		--j;
		while ( ++j < lgt ) {
		    iptn <<= B2; lptn <<= B4; 
		    switch ( seq.charAt(j) ) { 
		    case 'A':             lptn |= B1; break; case 'C': iptn |= B1; lptn |= B2; break; 
		    case 'G': iptn |= B2; lptn |= B4; break; case 'T': iptn |= B3; lptn |= B8; break; 
		    }
		    bsrkmer.set((iptn &= imsk)); lts.add( new Long((lptn &= lmsk)) );
		}
	    }
	    lrkmer = new long[(lrkmerLgt=(o=lts.size()))]; while ( --o >= 0 ) lrkmer[o] = lts.pollLast().longValue(); 
	    System.out.println(lrkmer.length + " rev k-mers (k=" + k + ")"); 
	    lts = null; contaminant = null; sb = null;

	    //################################
	    //## trimming reads (fwd & rev) ##
	    //################################
	    fin = new BufferedReader(new FileReader(finfile)); fout = new BufferedWriter(new FileWriter(foutfile));
	    rin = new BufferedReader(new FileReader(rinfile)); rout = new BufferedWriter(new FileWriter(routfile));
	    sout = new BufferedWriter(new FileWriter(soutfile));
	    pcpt = 0; ftcpt = 0; rtcpt = 0; rtcpt = 0; rrcpt = 0;
	    cpt = B0; score = new byte[51];
	    while ( true ) {
		//## reading fastq ##
		try { line = fin.readLine().trim(); } catch ( NullPointerException e ) { fin.close(); rin.close(); break; }
		switch ( ++cpt ) { case B1: fid1 = line; break; case B2: fseq = line.toUpperCase(); break; case B3: fid2 = line; break; case B4: fqsc = line; break; }
		try { line = rin.readLine().trim(); } catch ( NullPointerException e ) { fin.close(); rin.close(); break; }
		switch ( cpt ) { case B1: rid1 = line; continue; case B2: rseq = line.toUpperCase(); continue; case B3: rid2 = line; continue; case B4: rqsc = line; cpt = B0; break; }
		
		//## displaying info ##
		if ( ++pcpt % 1000000 == 0 ) 
		    System.out.println("[" + String.format(Locale.US, "%02d", new Long((cur=(System.currentTimeMillis()-t)/1000)/60))
				       + ":" + (line="0"+(cur%60)).substring(line.length()-2) + "]" 
				       + String.format(Locale.US, "%,12d", new Integer(pcpt)) + " read pairs processed:" 
				       + String.format(Locale.US, "%,12d", new Integer(ftcpt+rtcpt)) + " trimmed" 
				       + " (fwd:" + String.format(Locale.US, "%,12d", new Integer(ftcpt))
				       + "  rev:" + String.format(Locale.US, "%,12d", new Integer(rtcpt)) + ")"
				       + String.format(Locale.US, "%,12d", new Integer(frcpt+rrcpt)) + " removed" 
				       + " (fwd:" + String.format(Locale.US, "%,12d", new Integer(frcpt))
				       + "  rev:" + String.format(Locale.US, "%,12d", new Integer(rrcpt)) + ")");
		
		//## matching k-mers (fwd) ##
		if ( (lgt=(short)fseq.length()) >= score.length ) { score = new byte[++lgt]; --lgt; } else Arrays.fill(score, B0);
		end = lgt; start = B_1; start0 = B_1; scount = 0; iptn = B0; i = B_1; j = B_1; b = B0; bq = B0; 
		if ( (fseq.indexOf('N') != -1) || (fseq.indexOf('n') != -1) ) {
		    //## ... with degenerate nucleotides ##
		    lptn = L0; nptn = L0; n = B0;
		    while ( ++j < k_1 ) { 
			iptn <<= B2; lptn <<= B4; nptn <<= B4; 
			switch ( fseq.charAt(j) ) { 
			case 'A': case 'a':             lptn |= B1; nptn |= B1;  break; 
			case 'C': case 'c': iptn |= B1; lptn |= B2; nptn |= B2;  break; 
			case 'G': case 'g': iptn |= B2; lptn |= B4; nptn |= B4;  break; 
			case 'T': case 't': iptn |= B3; lptn |= B8; nptn |= B8;  break; 
			default:  n = k;                            nptn |= B15; 
			} 
			--n;
		    }
		    while ( j < lgt ) {
			b = (score[++i] += b); bq = ( fqsc.charAt(i) < qmin ) ? B1 : B0;
			iptn <<= B2; lptn <<= B4; nptn <<= B4; 
			switch ( fseq.charAt(j) ) { 
			case 'A': case 'a':             lptn |= B1; nptn |= B1;  break; 
			case 'C': case 'c': iptn |= B1; lptn |= B2; nptn |= B2;  break; 
			case 'G': case 'g': iptn |= B2; lptn |= B4; nptn |= B4;  break; 
			case 'T': case 't': iptn |= B3; lptn |= B8; nptn |= B8;  break; 
			default:  n = k;                            nptn |= B15; 
			} 
			switch ( (((--n) < B0) ? B0 : B1) ) {
			case B0:
			    if ( bsfkmer.get((iptn &= imsk)) ) { 
				b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
				scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
				start = ( (j >= start0) && (scount+k2 >= j) ) 
				    ? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				    : ++start; --start;
				continue; 
			    }
			    break;
			default:
			    id = Arrays.binarySearch(lfkmer, (lptn &= lmsk)); id = -(++id); --id;
			    o = Arrays.binarySearch(lfkmer, (nptn &= lmsk)); o = -(++o); 
			    while ( (++id < o) && ((lkm=lfkmer[id]) != (lkm |= (lpt=lptn))) ) {}
			    if ( id != o ) { 
				b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
				scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
				start = ( (j >= start0) && (scount+k2 >= j) ) 
				    ? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				    : ++start; --start;
				continue; 
			    }
			    break;
			}
			b = (score[i]=((b < B0) ? B0 : b)); ++j; 
			if ( (score[i] += bq) == B0 ) continue;
			scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
			start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		    }
		}
		else { 
		    //## ... with no degenerate nucleotide ##
		    while ( ++j < k_1 ) { 
			iptn <<= B2; 
			switch ( fseq.charAt(j) ) { case 'C': case 'c': iptn |= B1; break; case 'G': case 'g': iptn |= B2; break; case 'T': case 't': iptn |= B3; break; } 
		    }
		    while ( j < lgt ) {
			b = (score[++i] += b); bq = ( fqsc.charAt(i) < qmin ) ? B1 : B0;
			iptn <<= B2; switch ( fseq.charAt(j) ) { case 'C': case 'c': iptn |= B1; break; case 'G': case 'g': iptn |= B2; break; case 'T': case 't': iptn |= B3; break; } 
			if ( bsfkmer.get((iptn &= imsk)) ) { 
			    b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
			    scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
			    start = ( (j >= start0) && (scount+k2 >= j) ) 
				? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				: ++start; --start;
			    continue; 
			}
			b = (score[i]=((b < B0) ? B0 : b)); ++j; 
			if ( (score[i] += bq) == B0 ) continue;
			scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
			start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		    }
		}
		while ( ++i < lgt ) { 
		    b = (score[i] += b); bq = ( fqsc.charAt(i) < qmin ) ? B1 : B0; b = (score[i]=((b < B0) ? B0 : b)); 
		    if ( (score[i] += bq) == B0 ) continue;
		    scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
		    start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		}
		start = ( (start0 > B0) && (start0 >= lgt-start) ) ? B_1 : start;

		//## looking for 3' trimming index ##
		if ( (start0 >= B0) && (lgt-start-B1 >= minLgt) ) {
		    end0 = end; end = lgt; 
		    j = end0; i = j; i -= k; ++i; 
		    ++end0; end0 *= B2;
		    scount = B2;
		    end = ( (end > j) && (end-j <= mismatch) && (end0-j >= lgt) && (scount >= lgt-j) ) ? j : end;
		    start0 /= B2; o = ( start0 > start ) ? start0 : start; 
		    while ( (x=(--j)) >= o ) {
			--i; 
			if ( (b=score[j]) == B0 ) continue;
			scount += B2;
			end = ( (i >= B0) && (b >= B2) && ((--b) == score[++x]) && (score[i] != B0) ) 
			      ? ( (end0-i >= lgt) && (scount+k2 >= lgt-i) ) 
			        ? (( end == lgt ) ? i : (( end-j <= mismatch ) ? i : end)) 
			        : end
			      : end ;
			end = ( (end > j) && (end-j <= mismatch) && (end0-j >= lgt) && (scount >= lgt-j) ) ? j : end;
		    }
		    end = ( start < end ) ? --end : start; ++end;
		}
		else end = lgt;
		// here the trimmed read is read.substring(++start, end)

		//## displaying results (fwd) ##
		if ( (verbose != B0) && ((start != B_1) || (end != lgt)) ) {
		    System.out.println(fid1); System.out.println(fseq);
		    o = B_1; while ( ++o < lgt ) System.out.print( (( (b=score[o]) == B0 ) ? " " : ( b < B10 ) ? b : ( b == B10 ) ? "0" : (""+(char)(54+b))) ); System.out.println("");
		    o = B_1; while ( ++o <= start ) System.out.print(">"); --o; while ( ++o < end ) System.out.print(" ");  --o; while ( ++o < lgt ) System.out.print("<"); 
		    System.out.println(""); System.out.println("");
		}

		// verifying composition of remaining incorrect nucleotides
		l = end; l -= start; --l;
		scount = l; if ( prop != B0 ) { i = start; while ( ++i < end ) scount -= ( score[i] == B1 ) ? B1 : B0; }

		//## updating read and q-score (fwd) ##
		if ( (l < minLgt) || (100*scount < prop*l) ) { ++frcpt; fseq = ""; fqsc = ""; }
		else if ( (start != B_1) || (end != lgt) ) { ++ftcpt; fseq = fseq.substring(++start,end); fqsc = fqsc.substring(start--,end); }
		
		//## matching k-mers (rev) ##
		if ( (lgt=(short)rseq.length()) >= score.length ) { score = new byte[++lgt]; --lgt; } else Arrays.fill(score, B0);
		end = lgt; start = B_1; start0 = B_1; scount = 0; iptn = B0; i = B_1; j = B_1; b = B0; bq = B0; 
		if ( (rseq.indexOf('N') != -1) || (rseq.indexOf('n') != -1) ) {
		    //## ... with degenerate nucleotides ##
		    lptn = L0; nptn = L0; n = B0;
		    while ( ++j < k_1 ) { 
			iptn <<= B2; lptn <<= B4; nptn <<= B4; 
			switch ( rseq.charAt(j) ) { 
			case 'A': case 'a':             lptn |= B1; nptn |= B1;  break; 
			case 'C': case 'c': iptn |= B1; lptn |= B2; nptn |= B2;  break; 
			case 'G': case 'g': iptn |= B2; lptn |= B4; nptn |= B4;  break; 
			case 'T': case 't': iptn |= B3; lptn |= B8; nptn |= B8;  break; 
			default:  n = k;                            nptn |= B15; 
			} 
			--n;
		    }
		    while ( j < lgt ) {
			b = (score[++i] += b); bq = ( rqsc.charAt(i) < qmin ) ? B1 : B0;
			iptn <<= B2; lptn <<= B4; nptn <<= B4; 
			switch ( rseq.charAt(j) ) { 
			case 'A': case 'a':             lptn |= B1; nptn |= B1;  break; 
			case 'C': case 'c': iptn |= B1; lptn |= B2; nptn |= B2;  break; 
			case 'G': case 'g': iptn |= B2; lptn |= B4; nptn |= B4;  break; 
			case 'T': case 't': iptn |= B3; lptn |= B8; nptn |= B8;  break; 
			default:  n = k;                            nptn |= B15; 
			} 
			switch ( (((--n) < B0) ? B0 : B1) ) {
			case B0:
			    if ( bsrkmer.get((iptn &= imsk)) ) { 
				b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
				scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
				start = ( (j >= start0) && (scount+k2 >= j) ) 
				    ? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				    : ++start; --start;
				continue; 
			    }
			    break;
			default:
			    id = Arrays.binarySearch(lrkmer, (lptn &= lmsk)); id = -(++id); --id;
			    o = Arrays.binarySearch(lrkmer, (nptn &= lmsk)); o = -(++o); 
			    while ( (++id < o) && ((lkm=lrkmer[id]) != (lkm |= (lpt=lptn))) ) {}
			    if ( id != o ) { 
				b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
				scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
				start = ( (j >= start0) && (scount+k2 >= j) ) 
				    ? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				    : ++start; --start;
				continue; 
			    }
			    break;
			}
			b = (score[i]=((b < B0) ? B0 : b)); ++j; 
			if ( (score[i] += bq) == B0 ) continue;
			scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
			start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		    }
		}
		else { 
		    //## ... with no degenerate nucleotide ##
		    while ( ++j < k_1 ) { 
			iptn <<= B2; 
			switch ( rseq.charAt(j) ) { case 'C': case 'c': iptn |= B1; break; case 'G': case 'g': iptn |= B2; break; case 'T': case 't': iptn |= B3; break; } 
		    }
		    while ( j < lgt ) {
			b = (score[++i] += b); bq = ( rqsc.charAt(i) < qmin ) ? B1 : B0;
			iptn <<= B2; switch ( rseq.charAt(j) ) { case 'C': case 'c': iptn |= B1; break; case 'G': case 'g': iptn |= B2; break; case 'T': case 't': iptn |= B3; break; } 
			if ( bsrkmer.get((iptn &= imsk)) ) { 
			    b = (score[i]=(++b)); score[i] += bq; --score[++j]; 
			    scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0; 
			    start = ( (j >= start0) && (scount+k2 >= j) ) 
				? (( start < B0 ) ? j : (( i-start <= mismatch ) ? j : ++start)) 
				: ++start; --start;
			    continue; 
			}
			b = (score[i]=((b < B0) ? B0 : b)); ++j; 
			if ( (score[i] += bq) == B0 ) continue;
			scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
			start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		    }
		}
		while ( ++i < lgt ) { 
		    b = (score[i] += b); bq = ( rqsc.charAt(i) < qmin ) ? B1 : B0; b = (score[i]=((b < B0) ? B0 : b)); 
		    if ( (score[i] += bq) == B0 ) continue;
		    scount += B2; start0 = ( start0 < B0 ) ? ((short)(B2*i)) : start0;
		    start = ( (start < i) && ((x=i)-start <= mismatch) && ((++x) >= start0) && (scount > x) ) ? i : start; end = i; 
		}
		start = ( (start0 > B0) && (start0 >= lgt-start) ) ? B_1 : start;

		//## looking for 3' trimming index ##
		if ( (start0 >= B0) && (lgt-start-B1 >= minLgt) ) {
		    end0 = end; end = lgt; 
		    j = end0; i = j; i -= k; ++i; 
		    ++end0; end0 *= B2;
		    scount = B2;
		    end = ( (end > j) && (end-j <= mismatch) && (end0-j >= lgt) && (scount >= lgt-j) ) ? j : end;
		    start0 /= B2; o = ( start0 > start ) ? start0 : start; 
		    while ( (x=(--j)) >= o ) {
			--i; 
			if ( (b=score[j]) == B0 ) continue;
			scount += B2;
			end = ( (i >= B0) && (b >= B2) && ((--b) == score[++x]) && (score[i] != B0) ) 
			      ? ( (end0-i >= lgt) && (scount+k2 >= lgt-i) ) 
			        ? (( end == lgt ) ? i : (( end-j <= mismatch ) ? i : end)) 
			        : end
			      : end ;
			end = ( (end > j) && (end-j <= mismatch) && (end0-j >= lgt) && (scount >= lgt-j) ) ? j : end;
		    }
		    end = ( start < end ) ? --end : start; ++end;
		}
		else end = lgt;
		// here the trimmed read is read.substring(++start, end)

		//## displaying results (rev) ##
		if ( (verbose != B0) && ((start != B_1) || (end != lgt)) ) {
		    System.out.println(rid1); System.out.println(rseq);
		    o = B_1; while ( ++o < lgt ) System.out.print( (( (b=score[o]) == B0 ) ? " " : ( b < B10 ) ? b : ( b == B10 ) ? "0" : (""+(char)(54+b))) ); System.out.println("");
		    o = B_1; while ( ++o <= start ) System.out.print(">"); --o; while ( ++o < end ) System.out.print(" ");  --o; while ( ++o < lgt ) System.out.print("<"); 
		    System.out.println(""); System.out.println("");
		}

		// verifying composition of remaining incorrect nucleotides
		l = end; l -= start; --l;
		scount = l; if ( prop != B0 ) { i = start; while ( ++i < end ) scount -= ( score[i] == B1 ) ? B1 : B0; }

		//## updating read and q-score (rev) ##
		if ( (l < minLgt) || (100*scount < prop*l) ) { ++rrcpt; rseq = ""; rqsc = ""; }
		else if ( (start != B_1) || (end != lgt) ) { ++rtcpt; rseq = rseq.substring(++start,end); rqsc = rqsc.substring(start--,end); }

		//## writing trimmed read ##
		if ( (fseq.length() >= minLgt) && (rseq.length() >= minLgt) ) {
		    fout.write(fid1); fout.newLine(); fout.write(fseq); fout.newLine(); fout.write(fid2); fout.newLine(); fout.write(fqsc); fout.newLine();
		    rout.write(rid1); rout.newLine(); rout.write(rseq); rout.newLine(); rout.write(rid2); rout.newLine(); rout.write(rqsc); rout.newLine();
		    continue;
		}
		if ( (fseq.length() >= minLgt) && (rseq.length() < minLgt) ) {
		    sout.write(fid1); sout.newLine(); sout.write(fseq); sout.newLine(); sout.write(fid2); sout.newLine(); sout.write(fqsc); sout.newLine();
		    continue;
		}
		if ( (fseq.length() < minLgt) && (rseq.length() >= minLgt) ) {
		    sout.write(rid1); sout.newLine(); sout.write(rseq); sout.newLine(); sout.write(rid2); sout.newLine(); sout.write(rqsc); sout.newLine();
		}
	    }	
	    fout.close(); rout.close(); sout.close(); 
	    System.out.println("[" + String.format(Locale.US, "%02d", new Long((cur=(System.currentTimeMillis()-t)/1000)/60))
			       + ":" + (line="0"+(cur%60)).substring(line.length()-2) + "]" 
			       + String.format(Locale.US, "%,12d", new Integer(pcpt)) + " read pairs processed:" 
			       + String.format(Locale.US, "%,12d", new Integer(ftcpt+rtcpt)) + " trimmed" 
			       + " (fwd:" + String.format(Locale.US, "%,12d", new Integer(ftcpt))
			       + "  rev:" + String.format(Locale.US, "%,12d", new Integer(rtcpt)) + ")"
			       + String.format(Locale.US, "%,12d", new Integer(frcpt+rrcpt)) + " removed" 
			       + " (fwd:" + String.format(Locale.US, "%,12d", new Integer(frcpt))
			       + "  rev:" + String.format(Locale.US, "%,12d", new Integer(rrcpt)) + ")");
	}



    }

}

