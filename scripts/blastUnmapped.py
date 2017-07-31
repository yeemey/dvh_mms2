#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 16:31:59 2016
python 3.5.2
BioPython 1.68

@author: ymseah
"""
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Blast import NCBIWWW, NCBIXML
import re, os, subprocess, time

class run_software:
    
    def get_sample_name(self, fastq_file):
        """
        Get sample names from fastq files that end in '1.fastq', '2.fastq', 
        or '_sing.fastq'.
        """
        if re.search('[1-2].fastq|_[1-2].fastq|_R[1-2].fastq|_sing.fastq', fastq_file):
            sample_name_end_index = re.search('[1-2].fastq|_[1-2].fastq|_R[1-2].fastq|_sing.fastq|__sing.fastq|_R_sing.fastq', fastq_file).start()
            if re.search('/', fastq_file):
                start_indices = [m.start() for m in re.finditer('/', fastq_file)]
                start_indices.sort()
                sample_name_start_index = start_indices[-1] + 1
            else:
                sample_name_start_index = 0
            sample_name = fastq_file[sample_name_start_index:sample_name_end_index]
            return sample_name
        else:
            print(fastq_file)
            print('Input file name must end in \'_1.fastq\', \'_2.fastq\', \'_R1.fastq\', or \'_R2.fastq\'.')
    
    def get_all_sample_names(self, sample_dir):
        sample_files = os.listdir(sample_dir)
        sample_files.sort()
        sample_names = []
        for sample in sample_files:
            sample_name = self.get_sample_name(sample)
            if sample_name != None:
                sample_names.append(sample_name)
        return set(sample_names)
    
    def run_pear(self, fastq_f, fastq_r, output_dir, sample_name):
        """
        PEAR usage: -f <fastq_f>, -r <fastq_r> -o <output_base_name>
        """
        args = ['/Users/ymseah/Applications/pear-0.9.10-bin-64/./pear-0.9.10-bin-64', '-f', fastq_f, '-r', fastq_r, '-o', output_dir+sample_name]
        rp = subprocess.Popen(args)
        #allow PEAR to be KeyboardInterrupted
        try:
            while rp.poll() is None:
                time.sleep(0.1)
        except KeyboardInterrupt:
            rp.kill()
            raise
    
    def batch_run_pear(self, sample_dir, output_dir):
        sample_files = os.listdir(sample_dir)
        sample_files.sort()
        print('Getting sample names from FASTQ file names...')
        all_names = self.get_all_sample_names(sample_dir)
        fastq_f = []
        fastq_r = []
        for name in all_names:
            if name + '1.fastq' in sample_files:
                fastq_f.append(name + '1.fastq')
            elif name + '_1.fastq' in sample_files:
                fastq_f.append(name + '_1.fastq')
            elif name + '_R1.fastq' in sample_files:
                fastq_f.append(name + '_R1.fastq')
        for name in all_names:
            if name + '2.fastq' in sample_files:
                fastq_r.append(name + '2.fastq')
            elif name + '_2.fastq' in sample_files:
                fastq_r.append(name + '_2.fastq')
            elif name + '_R2.fastq' in sample_files:
                fastq_r.append(name + '_R2.fastq')
        
        fastq_f.sort()
        fastq_r.sort()
        samples = list(all_names)
        samples.sort()
        
        #check lengths of fastq_f and fastq_r
        if len(fastq_f) != len(fastq_r):
            print('Missing either forward/reverse read file(s)...')
        else:
            print('Equal number of forward and reverse read files.')
            #check if all samples have read files
            if len(samples) != len(fastq_f):
                print('Missing read files for samples...')
            else:
                sample_counter = 0
                while sample_counter < len(samples):
                    if re.search(samples[sample_counter], fastq_f[sample_counter]) and re.search(samples[sample_counter], fastq_r[sample_counter]):
                        print('Running PEAR on ' + samples[sample_counter])
                        #self.run_pear(sample_dir+fastq_f[sample_counter], sample_dir+fastq_r[sample_counter], output_dir, samples[sample_counter])
                    sample_counter += 1    
        
    def get_all_pear_samples(self, pear_output_dir):
        pear_files = os.listdir(pear_output_dir)
        all_pear_samples = []
        for file in pear_files:
            filenamesplit = file.split('.')
            if filenamesplit[1] == 'assembled' or filenamesplit[1] == 'unassembled' or filenamesplit[1] == 'discarded':
                all_pear_samples.append(filenamesplit[0])
            else:
                all_pear_samples.append(filenamesplit[0] + '.' + filenamesplit[1])
        return set(all_pear_samples)
        
    def run_breseq(self, input_dir, sample_name, output_dir, polymorphism_min = '5 ', ref_dir = '/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/data/', ref_genome1 = 'dv.gbk', ref_genome2 = 'mp.gbk', ref_genome3 = 'megaplasma.gbk'):
        """
        breseq usage: breseq -p -o . -r <ref_genome> --polymorphism-minimum-coverage-each-strand 5 <input_read_file>
        """
        breseq_dir = output_dir + sample_name + '_breseq'
        os.mkdir(breseq_dir)
        ref1 = ref_dir + ref_genome1
        ref2 = ref_dir + ref_genome2
        ref3 = ref_dir + ref_genome3
        reads1 = input_dir + sample_name + '.assembled.fastq'
        reads2 = input_dir + sample_name + '.unassembled.forward.fastq'
        reads3 = input_dir + sample_name + '.unassembled.reverse.fastq'
        
        args = ['breseq', '-p', '-o', breseq_dir, '-r', ref1, '-r', ref2, '-r',
                ref3, '--polymorphism-minimum-coverage-each-strand', polymorphism_min, 
                reads1, reads2, reads3]
        
        rb = subprocess.Popen(args)
        #allow breseq to be KeyboardInterrupted
        try:
            while rb.poll() is None:
                time.sleep(0.1)
        except KeyboardInterrupt:
            rb.kill()
            raise

    def batch_run_breseq(self, pear_results_dir, breseq_output_dir):
        all_samples = list(self.get_all_pear_samples(pear_results_dir))
        for sample in all_samples:
            self.run_breseq(pear_results_dir, sample, breseq_output_dir)
    
    def run_gdtools_compare(self, *samples, breseq_dir = '/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/breseq_results/', 
                            ref_dir = '/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/data/', 
                            ref_genome1 = 'dv.gbk', ref_genome2 = 'mp.gbk', ref_genome3 = 'megaplasma.gbk'):
        """
        Usage: gdtools COMPARE [-o annotated.html] -r reference.gbk input.1.gd [input.2.gd ... ]
        """
        ref1 = ref_dir + ref_genome1
        ref2 = ref_dir + ref_genome2
        ref3 = ref_dir + ref_genome3
        ancestor = breseq_dir + 'sic_Ancestor_breseq/output/0.gd'
        evo_line = samples[0]
        evo_line = evo_line[4:]
        evo_line_end_index = re.search('\.|-|_', evo_line).start()
        evo_line = evo_line[:evo_line_end_index]
        output_file = breseq_dir + 'compare/' + evo_line + '.html'
        args = ['gdtools', 'COMPARE', '-o', output_file, '-r', ref1, '-r', ref2,'-r', ref3, ancestor]

        sample_counter = 0
        while sample_counter < len(samples):
            gd_file = 'output.gd'
            gd_path = breseq_dir + samples[sample_counter] + '_breseq/output/'            
            if re.search('-13', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-100.gd'
                os.rename(gd_path + gd_file, new_gd_filepath)
            elif re.search('-14', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-100.gd'
                os.rename(gd_path + gd_file, new_gd_filepath)
            elif re.search('-15', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-100.gd'
                os.rename(gd_path + gd_file, new_gd_filepath) 
            elif re.search('-16', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-100.gd'
                os.rename(gd_path + gd_file, new_gd_filepath) 
            elif re.search('-18', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-100.gd'
                os.rename(gd_path + gd_file, new_gd_filepath)
            elif re.search('.45', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-300.gd'
                os.rename(gd_path + gd_file, new_gd_filepath)
            elif re.search('-76', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-500.gd'
                os.rename(gd_path + gd_file, new_gd_filepath)
            elif re.search('.118', samples[sample_counter]):
                new_gd_filepath = gd_path + evo_line + '-780.gd'
                os.rename(gd_path + gd_file, new_gd_filepath)
            else:
                new_gd_filepath = gd_path + evo_line + '-1000.gd'
                os.rename(gd_path + gd_file, new_gd_filepath)
            args.append(new_gd_filepath)
            sample_counter +=1
        
        rgc = subprocess.Popen(args)
        #allow gdtools to be KeyboardInterrupted
        try:
            while rgc.poll() is None:
                time.sleep(0.1)
        except KeyboardInterrupt:
            rgc.kill()
            raise

    def biopy_fastq_to_fasta(self, fastq_file, fasta_file):
        """Uses BioPython SeqIO module to convert fastq to fasta file"""
        SeqIO.convert(fastq_file, "fastq", fasta_file, "fasta")
        
    def batch_biopy_fastq_to_fasta(self, input_dir = '/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/breseq_results/'):
        fastq_filepath = input_dir #+ sample_name + '_breseq/data/'
        files = os.listdir(fastq_filepath)
        fastq_filenames = []
        fasta_filenames = []
        for each in files:
            if re.search('.unmatched.fastq', each):
                fastq_filenames.append(each)
        for fastq in fastq_filenames:
            self.biopy_fastq_to_fasta(fastq_filepath + fastq, fastq_filepath + fastq[:-5] + 'fasta')
            fasta_filenames.append(fastq[:-5] + 'fasta')
        return fasta_filenames
    
    def run_qblast_fasta(self, fasta_file):
        """Iterate through fasta file to individually BLAST each read"""
        for read in SeqIO.parse(fasta_file, "fasta"):
            result_handle = NCBIWWW.qblast("blastn", "nt", read.format("fasta"), hitlist_size=3)
            with open("my_blast.xml", "a") as blast_results:
                blast_results.write(result_handle.read())

    def run_blastn_remote(self, input_dir, fasta_file, output_dir):
        """Uses BLAST+ command line to search NCBI servers
        blastn usage: ./blastn -db nt -query <input_file> -out <output_file> -remote
        """
        blast_out = output_dir + fasta_file + 'BLASTout.txt'
        args = ['/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/bin/ncbi-blast-2.5.0+/bin/./blastn',
                '-db', 'nt', '-query', input_dir + fasta_file, '-out', blast_out, '-remote']
        rbn = subprocess.Popen(args)
        #allow blastn to be KeyboardInterrupted
        try:
            while rbn.poll() is None:
                time.sleep(0.1)
        except KeyboardInterrupt:
            rbn.kill()
            raise
    
    def batch_run_blastn_remote(self, list_of_files, input_dir = '/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/blast_results/input', 
                                blast_out_dir = '/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/blast_results/'):
        for fasta in fasta_filenames:
            self.run_blastn_remote(input_dir, fasta, blast_out_dir)

class parse_results:
    
    def parse_blast(self, blast_out):
        """Argument: blastn output text file name
        1. Records number of queries
        2. Records top (first) blast hit to D. vulgaris, M. maripaludis, & others.
        3. Records no hits found.
        4. Output: 2 files 
            a. Number of D.v., M.m., and no hits, with list of results.
            b. Number of other hits, with list of results.
            """
        with open(blast_out) as blast_results:
            file_lines = blast_results.readlines()
        num_queries = 0
        desulfo_v_hits = []
        methano_m_hits = []
        no_hits_counter = 0
        other_hits = []
        line_counter = 0
        while line_counter < len(file_lines):
            if re.match("Query=", file_lines[line_counter]):
                num_queries += 1
            elif re.search("Sequences producing significant alignments:", file_lines[line_counter]):
                top_hit_index = line_counter + 2
                if re.search("Desulfovibrio vulgaris", file_lines[top_hit_index]):
                    desulfo_v_hits.append(file_lines[top_hit_index])
                elif re.search("Methanococcus maripaludis", file_lines[top_hit_index]):
                    methano_m_hits.append(file_lines[top_hit_index])
                else:
                    other_hits.append(file_lines[top_hit_index])                
            elif re.search("No hits found", file_lines[line_counter]):
                no_hits_counter += 1
            line_counter += 1
        print("\nDesulfovibrio vulgaris hits / Total queries: " + str(len(desulfo_v_hits)) + "/" + str(num_queries) 
              + "\nMethanococcus maripaludis hits / Total queries: " + str(len(methano_m_hits)) + "/" + str(num_queries) 
              + "\nNo hits / Total queries: " + str(no_hits_counter) + "/" + str(num_queries) 
              + "\nOther hits / Total queries: " + str(len(other_hits)) + "/" + str(num_queries) + "\n")
        output_filename_string = blast_out + "_parsed.txt"
        output_other_hits_filename_string = blast_out + "_parsed_others.txt"
        with open(output_filename_string, "w") as outfile:
            outfile.write("Desulfovibrio vulgaris hits / Total queries: " + str(len(desulfo_v_hits)) + "/" + str(num_queries) 
                  + "\nMethanococcus maripaludis hits / Total queries: " + str(len(methano_m_hits)) + "/" + str(num_queries)
                  + "\nNo hits / Total queries: " + str(no_hits_counter) + "/" + str(num_queries) + "\n\n")
            for dv in desulfo_v_hits:
                outfile.write(dv)
            for mm in methano_m_hits:
                outfile.write(mm)
        with open(output_other_hits_filename_string, "w") as outfile2:
            outfile2.write("Other hits / Total queries: " + str(len(other_hits)) + "/" + str(num_queries) + "\n\n")
            for hits in other_hits:
                outfile2.write(hits)

#biopy_fastq_to_fasta("sic_M1_pear.assembled.unmatched.fastq","sic_M1_pear.assembled.unmatched.fasta")
#parse_blast("sic_M1_pear.assembled.unmatchedBLAST.out") 
#run_qblast_fasta("sic_HA3.45_1.unmatched.fasta")