#!/usr/bin/env python
# python /diskmnt/Projects/ccRCC_longread/longWGS/t2t/deepsomatic/get.supReads.py -b /diskmnt/Projects/ccRCC_longread/longWGS/t2t/aligned_bam/HT282_Tumor_phased.bam -v /diskmnt/Projects/ccRCC_longread/longWGS/t2t/deepsomatic/passedVCF/HT282_Tumor.deepsomatic.vcf -o /diskmnt/Projects/ccRCC_longread/longWGS/t2t/deepsomatic/test -q 20 -r /diskmnt/Projects/ccRCC_longread/00_analysis_Yuwei/chm13v2.0.fa
import argparse
import pysam
import os

def extract_supporting_reads_ids(bam_file, vcf_file, output_dir, min_base_quality=20):
    """
    Extract IDs and HP tags of reads supporting SNVs or indels from a BAM file based on variants in a VCF file.
    
    Parameters:
    -----------
    bam_file : str
        Path to the input BAM file
    vcf_file : str
        Path to the input VCF file containing variants
    output_dir : str
        Directory to save read IDs and HP tags
    min_base_quality : int
        Minimum base quality to consider a read as supporting a variant
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Open BAM file
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except (ValueError, OSError) as e:
        print(f"Error opening BAM file: {e}")
        return
    
    # Check if reference genome file is available (for better insertion detection)
    reference = None
    if hasattr(bam, 'reference_filename') and bam.reference_filename:
        try:
            reference = pysam.FastaFile(bam.reference_filename)
            print(f"Using reference genome from: {bam.reference_filename}")
        except (ValueError, OSError) as e:
            print(f"Warning: Could not open reference genome file: {e}")
    
    # Open VCF file
    try:
        vcf = pysam.VariantFile(vcf_file)
    except (ValueError, OSError) as e:
        print(f"Error opening VCF file: {e}")
        bam.close()
        if reference:
            reference.close()
        return
    
    # Process each variant in the VCF file
    for variant in vcf:
        # Format variant identifier
        variant_id = f"{variant.chrom}_{variant.pos}_{variant.ref}_{','.join(variant.alts)}"
        
        # Determine if variant is SNV or indel
        is_deletion = any(len(alt) < len(variant.ref) for alt in variant.alts)
        is_insertion = any(len(alt) > len(variant.ref) for alt in variant.alts)
        variant_type = "SNV" if (not is_deletion and not is_insertion) else "INDEL"
        
        # Output file for supporting reads
        output_file = os.path.join(output_dir, f"{variant_id}_{variant_type}_reads.txt")
        
        # Calculate the appropriate region to fetch reads
        fetch_start = max(0, variant.pos - 10)  # Start a bit before variant to catch reads
        fetch_end = variant.pos + max(len(variant.ref), max(len(alt) for alt in variant.alts)) + 10
        
        # Fetch reads overlapping the variant position
        try:
            reads = bam.fetch(variant.chrom, fetch_start, fetch_end)
        except ValueError as e:
            print(f"Error fetching reads for {variant_id}: {e}")
            continue
        
        # If we have a reference genome and this is an insertion, get downstream sequence
        downstream_seq = ""
        if reference and is_insertion:
            try:
                # Get sequence downstream of insertion point to compare with insertion
                # VCF position is 1-based, but pysam is 0-based
                downstream_seq = reference.fetch(variant.chrom, variant.pos, variant.pos + 10).upper()
            except Exception as e:
                print(f"Warning: Could not fetch reference sequence for {variant_id}: {e}")
        
        supporting_reads = []
        
        # Process each read
        for read in reads:
            if read.is_unmapped or read.is_duplicate or read.is_qcfail or read.is_secondary:
                continue
                
            # For SNVs
            if variant_type == "SNV":
                # Get the aligned pairs to find the position in the read corresponding to the variant
                read_positions = read.get_aligned_pairs(matches_only=True)
                for read_pos, ref_pos in read_positions:
                    # Check if this position corresponds to the variant position
                    if ref_pos == variant.pos - 1:  # 0-based indexing in BAM
                        # Get the base at this position in the read
                        read_base = read.query_sequence[read_pos]
                        # Get the quality score for this base
                        base_qual = read.query_qualities[read_pos]
                        
                        # Check if the read base matches one of the alternate alleles
                        if read_base in variant.alts and base_qual >= min_base_quality:
                            # Get HP tag if present, otherwise set to "None"
                            hp_tag = read.get_tag("HP") if read.has_tag("HP") else "None"
                            supporting_reads.append((read.query_name, hp_tag))
                        break
            
            # For deletions
            elif is_deletion:
                # For deletions, check if read has deletion at the correct position
                read_start = read.reference_start  # 0-based
                
                # Skip reads that don't overlap the deletion
                if read_start > variant.pos - 1 or (read.reference_end is not None and read.reference_end < variant.pos):
                    continue
                
                # Parse CIGAR string to find deletions
                ref_pos = read_start
                for op, length in read.cigartuples:
                    # For deletions (operation code 2)
                    if op == 2:
                        # Calculate the start and end of this deletion in reference coordinates
                        del_start = ref_pos
                        del_end = ref_pos + length
                        
                        # For VCF format deletion (ATTC → A), the actual deletion is (TTC) starting at position pos+1
                        # Check if this deletion overlaps our variant's deletion
                        # The anchor base is at variant.pos-1 (0-based), so deletion starts at variant.pos
                        if del_start <= variant.pos and del_end >= variant.pos + len(variant.ref) - len(variant.alts[0]):
                            # Get HP tag if present, otherwise set to "None"
                            hp_tag = read.get_tag("HP") if read.has_tag("HP") else "None"
                            supporting_reads.append((read.query_name, hp_tag))
                            break
                    
                    # Update reference position for operations that consume reference
                    if op in [0, 2, 3, 7, 8]:  # M, D, N, =, X operations
                        ref_pos += length
            
            # For insertions
            elif is_insertion:
                # Handle cases like ref=A, alt=AT where the T might be the same as the base after A in reference
                read_start = read.reference_start  # 0-based
                
                # First, check if this read spans the insertion point
                if read_start > variant.pos - 1 or (read.reference_end is not None and read.reference_end < variant.pos):
                    continue
                
                # We need to analyze both the CIGAR string and the actual sequence
                ref_pos = read_start
                read_pos = 0
                
                # Create a list of (ref_pos, read_pos, op, length) entries for detailed analysis
                alignment = []
                for op, length in read.cigartuples:
                    # Skip clipped bases in read positioning
                    if op == 4:  # soft clip
                        read_pos += length
                        continue
                    if op == 5:  # hard clip
                        continue
                        
                    if op == 0 or op == 7 or op == 8:  # Match, Equal, or Mismatch
                        alignment.append((ref_pos, read_pos, op, length))
                        ref_pos += length
                        read_pos += length
                    elif op == 1:  # Insertion
                        alignment.append((ref_pos, read_pos, op, length))
                        read_pos += length
                    elif op == 2 or op == 3:  # Deletion or Skip
                        alignment.append((ref_pos, read_pos, op, length))
                        ref_pos += length
                
                # Analyze the alignment to find insertions
                supports_insertion = False
                for i, (r_pos, q_pos, op, length) in enumerate(alignment):
                    # Check for direct insertion operations at the variant position
                    if r_pos == variant.pos - 1 and op == 1:
                        # We have an insertion right at the variant position
                        inserted_bases = read.query_sequence[q_pos:q_pos+length]
                        
                        # For each alt allele, check if our insertion matches
                        for alt in variant.alts:
                            # The inserted sequence is alt[len(variant.ref):]
                            inserted_seq = alt[len(variant.ref):]
                            if inserted_bases == inserted_seq:
                                supports_insertion = True
                                break
                    
                    # Check for potential hidden insertions where the aligner might have matched
                    # For example, ref=A, alt=AT, downstream=T might be aligned as "match" instead of "insertion"
                    if not supports_insertion and i < len(alignment) - 1:
                        curr_pos, curr_q_pos, curr_op, curr_len = alignment[i]
                        next_pos, next_q_pos, next_op, next_len = alignment[i+1]
                        
                        # If we're at the variant position with a match
                        if curr_pos <= variant.pos - 1 < curr_pos + curr_len and curr_op in [0, 7, 8]:
                            # Position in the read sequence for the variant position
                            var_read_offset = curr_q_pos + (variant.pos - 1 - curr_pos)
                            
                            # For each alt allele
                            for alt in variant.alts:
                                inserted_seq = alt[len(variant.ref):]
                                
                                # Check if the read sequence contains the insertion
                                # We need to look for the entire expected sequence: ref base + insertion
                                expected_seq = variant.ref + inserted_seq
                                read_seq_to_check = read.query_sequence[var_read_offset:var_read_offset+len(expected_seq)]
                                
                                if read_seq_to_check == expected_seq:
                                    supports_insertion = True
                                    break
                
                if supports_insertion:
                    # Get HP tag if present, otherwise set to "None"
                    hp_tag = read.get_tag("HP") if read.has_tag("HP") else "None"
                    supporting_reads.append((read.query_name, hp_tag))
        
        # Write supporting read IDs and HP tags to file
        if supporting_reads:
            with open(output_file, "w") as out:
                out.write("Read_ID\tHP_Tag\n")
                for read_name, hp_tag in supporting_reads:
                    out.write(f"{read_name}\t{hp_tag}\n")
            
            # print(f"Found {len(supporting_reads)} supporting reads for {variant_id} ({variant_type})")
        # else:
            # print(f"No supporting reads found for {variant_id} ({variant_type})")
    
    # Clean up
    bam.close()
    vcf.close()
    if reference:
        reference.close()
    print(f"Extraction complete. Results saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Extract read IDs and HP tags supporting SNVs or indels")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-q", "--min-quality", type=int, default=20, help="Minimum base quality (default: 20)")
    parser.add_argument("-r", "--reference", help="Reference genome FASTA file (optional, improves insertion detection)")
    
    args = parser.parse_args()
    
    extract_supporting_reads_ids(args.bam, args.vcf, args.output, args.min_quality)

if __name__ == "__main__":
    main()
