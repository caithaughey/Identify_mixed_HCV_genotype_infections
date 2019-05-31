import os
import sys
import subprocess
import statistics
from Bio import SeqIO
from Bio import AlignIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from collections import Counter

#####################################################################################################
###
### Author: Caitlin Vigetun Haughey
###
### Usage: Run command:
###         python parse_blast_results.py 'folder/directory containing the reads file' 'threshold value for consensus sequence'
###
### Description: Checks the reads for patients HCV sample sequenced with PacBio (or possibly other long read technology)
###              to assess if the sample contains a mixed infection (an infection with multiple HCV genotypes).
###              Reads are then filtered based on length and number of gaps. The reads are aligned and a consensus
###              sequence is created for the sample. This consensus sequence is in fasta format and can be used
###              in the Geno2pheno database (https://hcv.geno2pheno.org) for RAV identification and genotyping.
###
#####################################################################################################


# Open the genotype reference sequence files (1a, 1b, 3a & 4a).
try:
    gt1a_reference_file = open("HCV_references/hcv-NS5A-gt1a.fasta")
    gt1b_reference_file = open("HCV_references/hcv-NS5A-gt1b.fasta")
    gt2b_reference_file = open("HCV_references/hcv-NS5A-gt2b.fasta")
    gt3a_reference_file = open("HCV_references/hcv-NS5A-gt3a.fasta")
    gt4a_reference_file = open("HCV_references/hcv-NS5A-gt4a.fasta")
except FileNotFoundError as fnf_error:
    print(fnf_error)
else:
    # Get the input file path
    query_file = sys.argv[1]
    query_folder = os.path.dirname(os.path.realpath(sys.argv[1]))
    query_folder = query_folder + "/"

    # To extract the file name, get the string between the directory slash ("/") and the file extension period (".").
    start = query_file.find('/', 0)
    end = query_file.find('.', start)
    file_name = query_file[start+1:end]

    if not os.path.isfile(query_file):
        print("\n")
        print('Input file does not exist, analysis terminated.')
        print("\n")
        sys.exit()
    else:
        try:
            threshold = float(sys.argv[2])
        except ValueError:
            print("\n")
            print("Input threshold value is not a number. Make sure to use a period decimal (.) and not a comma (,).")
            print("\n")
        else:
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            print("\n")
            print("Running analysis for sample " + file_name + "...")
            print("\n")


            ## ############################################################################
            ##
            ## If reads are in fastq format, convert and create a fasta file with the reads.
            ##
            ## ############################################################################

            # If the original fastq input file with the reads has not been converted into fasta, run SeqIO.convert to get a fasta version of the input.
            if not any(fname.endswith('.fasta') for fname in os.listdir(query_folder)):
                try:
                    open(query_folder + file_name + ".fastq")
                except FileNotFoundError as fnf_error:
                    print(fnf_error)
                    print("\n")
                    print("Check that the reads file exists in the given directory.")
                    print("\n")
                    sys.exit()
                else:
                    SeqIO.convert((query_folder + file_name + ".fastq"),'fastq',(query_folder + file_name + ".fasta"),'fasta')
                    print("Reads file successfully converted from fastq to fasta format!")
                    print("\n")

            # Save the relative path of the query file + the fasta file name.
            fasta_file_path = query_folder + file_name + ".fasta"

            try:
                os.makedirs(query_folder + "mixed_infection_files")
            except FileExistsError:
                # directory already exists
                pass
            try:
                os.makedirs(query_folder + "consensus_sequence_files")
            except FileExistsError:
                # directory already exists
                pass


            ## ###################################################################################
            ##
            ## Run reads through a local BLAST, against the NS5A genotype reference sequence files.
            ## Each read is matched against closest reference.
            ##
            ## ###################################################################################

            # If the input fasta reads file has not yet been BLASTed against the reference genotype sequences.
            if not any(fname.endswith('.xml') for fname in os.listdir(query_folder + "mixed_infection_files/")):
                # Save the blastn command with the input and output files, and with parameters for output file format, gap penalties, etc.
                command = ["blastn", "-db", "NS5A_genotypes", "-query", fasta_file_path, "-outfmt", "5", "-out", query_folder + "mixed_infection_files/" + file_name + ".xml", "-task", "blastn", "-reward", "1", "-dust", "yes", "-penalty", "-3", "-gapopen", "5", "-gapextend", "2", "-max_hsps", "1", "-max_target_seqs", "1"]
                # Pass the command as a subprocess to run through bash.
                try:
                    process = subprocess.run(command, stdout=subprocess.PIPE, check=True)
                    subprocess.CompletedProcess(args=command, returncode=0)
                except OSError:
                    print("ERROR: Blastn not found, make sure you have the BLAST+ software from ncbi installed.")
                    print("\n")
                    print("How to download, install and use BLAST+: https://www.ncbi.nlm.nih.gov/books/NBK279690/")
                    print("\n")
                    sys.exit()
                except subprocess.CalledProcessError as error:
                    print(error)
                    print("\n")
                    print('The nucleotide database "NS5A_genotypes" cannot be found in the blastdb. Make sure to create the NS5A_genotypes database, instructions can be found at: https://www.ncbi.nlm.nih.gov/books/NBK279688/.')
                    print("\n")
                    os.remove(fasta_file_path[:-6] + ".xml")
                    sys.exit()
                else:
                    print("Read files have been blasted against reference genotype sequences!")
                    print("\n")

            # Get the resulting .xml file from the blast and parse this file with NCBIXML.parse.
            file_list = [f for f in os.listdir(query_folder + "mixed_infection_files/") if f.endswith('.xml')]
            xml_file = file_list[0]
            result_handle = open(query_folder + "mixed_infection_files/" + xml_file)
            blast_records = NCBIXML.parse(result_handle)

            # Count the total number of reads in the input file.
            total_seq_count = 0
            with open(fasta_file_path) as fasta_file:
                for line in fasta_file:
                    if line.startswith(">"):
                        total_seq_count += 1
            fasta_file.close()


            ## ###################################################################################
            ##
            ## Sort the BLAST results for the reads based on genotype.
            ##
            ## ###################################################################################

            # Initialize empty counts for each genotype
            gt1a_seq_count = 0
            gt1b_seq_count = 0
            gt2b_seq_count = 0
            gt3a_seq_count = 0
            gt4a_seq_count = 0

            # Check if the reads have been sorted into genotype based total reads files, if not then sort the reads.
            if not any(fname.startswith(file_name + '_total_reads') for fname in os.listdir(query_folder + "mixed_infection_files/")):

                # Create empty dictionary and hits list.
                q_dict = SeqIO.index(fasta_file_path, "fasta")
                hits = []
                low_e_value_seq_count = 0

                # Set e-value threshold (any blast hit with an e-value above the threshold is ignored).
                E_VALUE_THRESH = 0.001
                # Loop through all blast records (for each read), and each alignment for each of those reads.
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        if alignment.hsps[0].expect <= E_VALUE_THRESH:
                            # Add the hit to the list.
                            hits.append(blast_record.query.split()[0])
                            # Check if the title of the alignment contains "GT1A", genotype 1a.
                            if "GT1A" in alignment.title:
                                gt1a_seq_count += 1
                                # Open (or create of doesn't exist yet) the total reads file for genotype 1a.
                                with open((query_folder + "mixed_infection_files/" + file_name + "_total_reads_gt1a.fasta"), "a") as gt1a:
                                    gt1a.write(">")
                                    # Write the read name to the file.
                                    gt1a.write(blast_record.query)
                                    gt1a.write("\n")
                                    # Get the read sequence from the alignment.
                                    hsp = alignment.hsps[0]
                                    seq = Seq(hsp.query)
                                    # If the read sequence does not start with "GAT", then the read is in reverse and is complemented.
                                    if not hsp.query.startswith("GAT"):
                                        seq = seq.reverse_complement()
                                    seq = str(seq)
                                    # Write the read sequence to the file.
                                    gt1a.write(seq)
                                    gt1a.write("\n")
                                gt1a.close()
                            # Same procedure as above, for genotype = gt1b.
                            if "GT1B" in alignment.title:
                                gt1b_seq_count += 1
                                with open((query_folder + "mixed_infection_files/" + file_name + "_total_reads_gt1b.fasta"), "a") as gt1b:
                                    gt1b.write(">")
                                    gt1b.write(blast_record.query)
                                    gt1b.write("\n")
                                    hsp = alignment.hsps[0]
                                    seq = Seq(hsp.query)
                                    if not hsp.query.startswith("GAT"):
                                        seq = seq.reverse_complement()
                                    seq = str(seq)
                                    gt1b.write(seq)
                                    gt1b.write("\n")
                                gt1b.close()
                            # Same procedure as above, for genotype = gt2b.
                            if "GT2B" in alignment.title:
                                gt2b_seq_count += 1
                                with open((query_folder + "mixed_infection_files/" + file_name + "_total_reads_gt2b.fasta"), "a") as gt2b:
                                    gt2b.write(">")
                                    gt2b.write(blast_record.query)
                                    gt2b.write("\n")
                                    hsp = alignment.hsps[0]
                                    seq = Seq(hsp.query)
                                    if not hsp.query.startswith("GAT"):
                                        seq = seq.reverse_complement()
                                    seq = str(seq)
                                    gt2b.write(seq)
                                    gt2b.write("\n")
                                gt2b.close()
                            # Same procedure as above, for genotype = gt3a.
                            if "GT3A" in alignment.title:
                                gt3a_seq_count += 1
                                with open((query_folder + "mixed_infection_files/" + file_name + "_total_reads_gt3a.fasta"), "a") as gt3a:
                                    gt3a.write(">")
                                    gt3a.write(blast_record.query)
                                    gt3a.write("\n")
                                    hsp = alignment.hsps[0]
                                    seq = Seq(hsp.query)
                                    if not hsp.query.startswith("GAT"):
                                        seq = seq.reverse_complement()
                                    seq = str(seq)
                                    gt3a.write(seq)
                                    gt3a.write("\n")
                                gt3a.close()
                            # Same procedure as above, for genotype = gt4a.
                            if "GT4A" in alignment.title:
                                gt4a_seq_count += 1
                                with open((query_folder + "mixed_infection_files/" + file_name + "_total_reads_gt4a.fasta"), "a") as gt4a:
                                    gt4a.write(">")
                                    gt4a.write(blast_record.query)
                                    gt4a.write("\n")
                                    hsp = alignment.hsps[0]
                                    seq = Seq(hsp.query)
                                    if not hsp.query.startswith("GAT"):
                                        seq = seq.reverse_complement()
                                    seq = str(seq)
                                    gt4a.write(seq)
                                    gt4a.write("\n")
                                gt4a.close()
                        else:
                            low_e_value_seq_count += 1
                            with open((query_folder + "mixed_infection_files/" + file_name + "_low_e_value.fasta"), "a") as low_file:
                                low_file.write(">")
                                # Write the read name to the file.
                                low_file.write(blast_record.query + "_" + alignment.title)
                                low_file.write("\n")
                                # Get the read sequence from the alignment.
                                hsp = alignment.hsps[0]
                                seq = Seq(hsp.query)
                                # If the read sequence does not start with "GAT", then the read is in reverse and is complemented.
                                if not hsp.query.startswith("GAT"):
                                    seq = seq.reverse_complement()
                                seq = str(seq)
                                # Write the read sequence to the file.
                                low_file.write(seq)
                                low_file.write("\n")
                            low_file.close()


                # Set number of misses as total dictionary of reads - reads that had atleast 1 blast hit.
                misses = set(q_dict.keys()) - set(hits)
                # If any missed reads, create "no match reads" fasta file and add these reads.
                if (len(misses) != 0):
                    with open((query_folder + "mixed_infection_files/" + file_name + "_no_match_reads.fasta"), "a+") as missed:
                        missed.write(str(misses))
                    missed.close()

                print("Blast results have been sorted by highest matching genotype reference!")
                print("\n")


            ## ###################################################################################
            ##
            ## Calculate the most common read length for the genotype sorted reads.
            ## This read length will be used to filter the reads by before creating a consensus sequence.
            ##
            ## ###################################################################################

            # Create empty list to hold the relative lengths of all reads for each specific total reads genotype file.
            len_list = []
            # Create empty list to hold the most common read length for each "total reads" per genotype file.
            all_seq_len_list = []

            # Create list of all "total reads" file names (in cases of mixed infections, there will be more than one "total reads" file).
            totals_file_list = [f for f in os.listdir(query_folder + "mixed_infection_files/") if 'total_reads' in f]
            # For each "total reads" file, extract the genotype from the file name.
            for i in range(len(totals_file_list)):
                file = totals_file_list[i]
                gt = file[-10:-6]
                # Open the "total reads" file.
                with open((query_folder + "mixed_infection_files/" + file), "r+") as total_file:
                    for line in total_file:
                        # For each line that is a read sequence line ( != > not a fasta header line).
                        if not line.startswith(">"):
                            # Append the length of the read sequence to the len_list.
                            len_list.append(len(line.strip()))
                            # Get the most common sequence length from the list.
                            c = Counter(len_list)
                            most_common_seq_len = float(c.most_common(1)[0][0])
                    # Save the genotype of these reads and the most common length for the specific "total reads" file as a tuple.
                    gt_len = (gt, most_common_seq_len)
                    # Append the genotype + read length tuple to the "all sequences length list".
                    all_seq_len_list.append(gt_len)
                total_file.close()
                len_list.clear()

            ## ###################################################################################
            ##
            ## Filter the reads based on whether or not they have the same length as the most common read length,
            ## number of gaps, and if the sequence starts with 'GAT' (the primer sequences first 3 positions).
            ##
            ## ###################################################################################

            # If the no filtered_reads file exists, create the file and filtered the reads.
            if not any(fname.startswith(file_name + '_filtered_reads') for fname in os.listdir(query_folder + "mixed_infection_files/")):
                # Loop through all the totals files (mixed infections will have more than one, since these are based on genotype).
                for i in range(len(totals_file_list)):
                    file = totals_file_list[i]
                    # Get the genotype of the totals file from the file name
                    gt = file[-10:-6]
                    # Open the totals file.
                    with open((query_folder + "mixed_infection_files/" + file), "r+") as total_file:
                        # Open/create a filtered reads file and name it after the sample id name and the genotype.
                        with open((query_folder + "consensus_sequence_files/" + file_name + "_filtered_reads_" + gt + ".fasta"), "a+") as filtered:
                            for line in total_file:
                                # If line starts with ">", save this as the header.
                                if line.startswith(">"):
                                    header_line = line
                                    # Get the sequence after the header using next().
                                    seq_line = next(total_file)
                                    # Count the number of gaps in the sequence ("-").
                                    nr_gaps = seq_line.count('-')
                                    # Create a list of all the genotypes in the most common sequence lengths list.
                                    gt_list = [x[0] for x in all_seq_len_list]
                                    # Get the index of the most common length for the current genotype (the genotype of the current totals file).
                                    j = gt_list.index(gt)
                                    # If the sequence has the same length as the most common read length for the genotype, has a maximum of 1 gap, and starts with "GAT".
                                    # The number of allowed gaps can be altered by the user by changing the value after "(nr_gaps < )".
                                    if ((len(seq_line.strip()) == all_seq_len_list[j][1]) and (nr_gaps < 2) and (seq_line.startswith("GAT"))):
                                        # If all criteria are met, add the header and the corresponding sequence to the filtered reads file.
                                        filtered.write(header_line)
                                        filtered.write(seq_line)
                        filtered.close()
                        print(gt + " reads have been filtered based on read length and number of gaps!")
                        print("\n")
                    total_file.close()


            ## ###################################################################################
            ##
            ## Calculate how many reads there are for each genotype,
            ## to evaluate whether or not the sample contains a mixed infection.
            ##
            ## ###################################################################################

            # If the "mixed_infection_results" file does not already exist, create the file and check if there is a mixed infection.
            if not any(fname.startswith(file_name + '_mixed_infection_results') for fname in os.listdir(query_folder + "mixed_infection_files/")):

                # Create an empty list to hold the frequencies of the genotypes (number of reads for genotype / number of total reads)
                freq_list = []

                # Calculate the frequencies for each genotype and append to the freq_list.
                gt1a_freq = gt1a_seq_count/total_seq_count
                freq_list.append(gt1a_freq)
                gt1b_freq = gt1b_seq_count/total_seq_count
                freq_list.append(gt1b_freq)
                gt2b_freq = gt2b_seq_count/total_seq_count
                freq_list.append(gt2b_freq)
                gt3a_freq = gt3a_seq_count/total_seq_count
                freq_list.append(gt3a_freq)
                gt4a_freq = gt4a_seq_count/total_seq_count
                freq_list.append(gt4a_freq)


                # Get the highest frequency from the list and the index of that frequency.
                highest_freq  = max(freq_list)
                highest_index = freq_list.index(highest_freq)
                # If the highest frequency is equal to or less than 0.85 (85%), this could indicate a mixed infection.
                if (highest_freq < 1):
                    # Create a copy of the freq_list and delete the highest frequency value.
                    freq_list2 = list(freq_list)
                    del freq_list2[-highest_index]
                    # Calculate the second highest frequency value in the copied freq list.
                    second_highest_freq  = max(freq_list2)
                    # If the second highest frequency is equal to or larger than 0.15 (15%), this is considered a mixed infection.
                    # These cut-off value (85% and 15%) can be altered by the user if they want the cut-off to be higher or lower.
                    if (second_highest_freq > 0):
                        mixed_infection = True
                        print("This sample contains a mixed genotype infection!")
                        print("\n")
                    else:
                        mixed_infection = False
                else:
                    mixed_infection = False

                # Open or create a text file with the mixed infection analysis results.
                with open((query_folder + "mixed_infection_files/" + file_name + "_mixed_infection_results.txt"), "w+") as totals:
                    # Add the count of the reads for each genotype, as well as the genotypes frequencies.
                    totals.write("Number of genotype 1a reads: ")
                    totals.write(str(gt1a_seq_count))
                    totals.write("\t")
                    totals.write("Frequency of genotype 1a reads: ")
                    totals.write(str(gt1a_freq))
                    totals.write("\n")
                    totals.write("Number of genotype 1b reads: ")
                    totals.write(str(gt1b_seq_count))
                    totals.write("\t")
                    totals.write("Frequency of genotype 1b reads: ")
                    totals.write(str(gt1b_freq))
                    totals.write("\n")
                    totals.write("Number of genotype 2b reads: ")
                    totals.write(str(gt2b_seq_count))
                    totals.write("\t")
                    totals.write("Frequency of genotype 2b reads: ")
                    totals.write(str(gt2b_freq))
                    totals.write("\n")
                    totals.write("Number of genotype 3a reads: ")
                    totals.write(str(gt3a_seq_count))
                    totals.write("\t")
                    totals.write("Frequency of genotype 3a reads: ")
                    totals.write(str(gt3a_freq))
                    totals.write("\n")
                    totals.write("Number of genotype 4a reads: ")
                    totals.write(str(gt4a_seq_count))
                    totals.write("\t")
                    totals.write("Frequency of genotype 4a reads: ")
                    totals.write(str(gt4a_freq))
                    totals.write("\n")
                    totals.write("\n")
                    totals.write("Total number of reads for this sample: ")
                    totals.write(str(total_seq_count))
                    totals.write("\n")
                    totals.write("\n")
                    totals.write("Number of reads with only low expected value match (e-value): ")
                    totals.write(str(low_e_value_seq_count))
                    totals.write("\n")
                    totals.write("\n")
                    totals.write("Number of reads with no NS5A genotype match: ")
                    totals.write(str(len(misses)))
                    totals.write("\n")

                    # If the mixed_infection = TRUE, the sample contains a mixed infection.
                    if (mixed_infection):
                        totals.write("\n")
                        totals.write("This is an infection with mixed genotypes.")
                    # Else the sample does not contain a mixed infection.
                    else:
                        totals.write("\n")
                        totals.write("This is not an infection with mixed genotypes.")
                totals.close()
                print("The mixed infection analysis has been completed!")
                print("\n")
                print("Total number of reads for this sample: ")
                print(str(total_seq_count))
                print("\n")
                print("Number of reads with only low expected value match (e-value): ")
                print(str(low_e_value_seq_count))
                print("\n")


            ## ###################################################################################
            ##
            ## Align the filtered read files using MUSCLE software.
            ##
            ##  * Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput
            ##  Nucleic Acids Res. 32(5):1792-1797.
            ##  * Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced time and space complexity
            ##  BMC Bioinformatics, (5) 113.
            ##
            ## ###################################################################################

            # If the "aligned_reads" file does not already exist, create the file and perform the alignment.
            if not any(fname.startswith(file_name + '_aligned_reads') for fname in os.listdir(query_folder + "consensus_sequence_files/")):
                # Create list of all the filtered reads files (mixed infections will have more than one).
                filtered_file_list = [f for f in os.listdir(query_folder + "consensus_sequence_files/") if 'filtered_reads' in f]
                # For each of the filtered reads files.
                for i in range(len(filtered_file_list)):
                    file = filtered_file_list[i]
                    # Get the genotype of the filtered file from the file name
                    gt = file[-10:-6]
                    # Save the path to the input file (the query folder + the filtered reads file)
                    input_path = (query_folder + "consensus_sequence_files/" + file)
                    # Save the path (and file name) of the resulting alignment file (the name of the file will contain the sample id + the genotype).
                    alignment_path = query_folder + "consensus_sequence_files/" + file_name + "_aligned_reads_" + gt + ".fasta"
                    # Check that the filtered reads file is not empty before running alignment.
                    if (os.path.getsize(input_path) > 0):
                        # Save the command to run muscle.

                        command = ["muscle3.8.31_i86darwin64", "-in", input_path, "-out", alignment_path]
                        # Pass the command as a subprocess to run through bash.
                        try:
                            process = subprocess.run(command, stdout=subprocess.PIPE, check=True)
                            subprocess.CompletedProcess(args=command, returncode=0)
                        except OSError:
                            print("ERROR: MUSCLE not found, make sure you have the software installed.")
                            print("\n")
                            print("How to download, install and use MUSCLE: https://www.drive5.com/muscle/")
                            print("\n")
                            sys.exit()
                        else:
                            print("Alignment using MUSCLE has been completed!")
                            print("\n")
                    # If input file is empty, skip this file.
                    else:
                        print("Input file empty, skipping file...")
                        print("\n")


            ## ###################################################################################
            ##
            ## Create a consensus sequence for the sample using the filtered reads.
            ## If a nucleotide position in the sequences differ, the threshold value chosen when the program was called is used.
            ## For example, if a position has A's for 87% of the reads and 13% G and a threshold of 0.85 was chosen,
            ## The consensus will be A ( 87 > 85). If the A's have 60% and G's 40%, the consensus will be the ambiguity
            ## character 'N' ( 60 < 85).
            ##
            ## ###################################################################################

            # If the "consensus.fasta" file does not already exist, create the file and create the consensus.
            if not any(("consensus" in fname) for fname in os.listdir(query_folder + "consensus_sequence_files/")):
                # Get the threshold value chosen by the user.
                threshold = float(sys.argv[2])
                # Create list of all the aligned reads files (mixed infections will have more than one).
                aligned_file_list = [f for f in os.listdir(query_folder + "consensus_sequence_files/") if 'aligned_reads' in f]
                # For each of the aligned reads files.
                for i in range(len(aligned_file_list)):
                    file = aligned_file_list[i]
                    # Get the genotype from the file name
                    gt = file[-10:-6]
                    # Pass the aligned reads file to the AlignIO module.
                    alignment = AlignIO.read(open(query_folder + "consensus_sequence_files/" + file_name + "_aligned_reads_" + gt + ".fasta"), 'fasta')
                    # Calculate summary info about the alignment.
                    summary_align = AlignInfo.SummaryInfo(alignment)
                    # Open or create a fasta file with the consensus sequence.
                    with open((query_folder + "consensus_sequence_files/" + file_name + "_consensus_" + gt + ".fasta"), "w+") as cons:
                        # Write a header to the sequence containing the sample id and the genotype.
                        cons.write(">Consensus_for_" + file_name + "_" + gt)
                        cons.write("\n")
                        # Create the consensus sequence using the function 'dumb_consensus', the threshold value and the ambiguity symbol 'N'.
                        # Write the consensus sequence to the consensus fasta file.
                        cons.write(str(summary_align.dumb_consensus(threshold, ambiguous="N")))
                    cons.close()
                    print("A consensus sequence has been created using the alignment file and the chosen threshold value!")
                    print("\n")
            # If a consensus sequence already exists, all steps of the program have already been run and all output files should exist.
            else:
                print("All files already exist for sample " + file_name + ", skipping...")
                print("\n")
