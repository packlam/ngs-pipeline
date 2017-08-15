import argparse
import getpass
import logging
import os
import shutil # for file copying
import subprocess
import sys
from time import strftime
import pandas as pd


# directories
ANALYSIS_DIR = '/HMDS_BIOINFORMATICS/HMDS/HMDS/MiSeq/Analysis/'
BACKUP_DIR = '/HMDS_BIOINFORMATICS/HMDS/HMDS/MiSeq/Backup/DNAseq/'
REFERENCES_DIR = '/HMDS_BIOINFORMATICS/references/'

# software
ABRA = '/usr/local/bfx/abra/abra2-2.07/abra2-2.07.jar'
ALAMUT = '/usr/local/bfx/alamut-batch/alamut-batch-standalone-1.7.0/alamut-batch'
BWA = '/usr/local/bfx/bwa/bwakit-0.7.12/bwa'
FASTQC = '/usr/local/bfx/fastqc/fastqc-0.11.5/fastqc'
PICARD = '/usr/local/bfx/picard/picard-2.9.4/picard.jar'
SAMTOOLS = '/usr/local/bfx/samtools/samtools-1.5/bin/samtools'
VARSCAN = '/usr/local/bfx/varscan/varscan-2.4.3/VarScan.v2.4.3.jar'

# python scripts
FILTERVCF = '/usr/local/bfx/python/ngs-pipeline/filtervcf.py'
GENERATE_RESULTS = '/usr/local/bfx/python/ngs-pipeline/generate_results.py'

# reference files
REF_GENOME = '/HMDS_BIOINFORMATICS/references/hg19/hg19.fa'
ASXL1_BED = '/HMDS_BIOINFORMATICS/references/myeloid/ASXL1.bed'


def get_args():
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = 'Runs the NGS pipeline',
        epilog = '''--run_id refers to the long name for the MiSeq run, e.g.
        \n\t170511_M02594_0240_000000000-B4YV9
        \n--expt_id refers to the short name for the run, e.g.
        \n\tMiSeq_myeloid_2017_19
        '''
    )
    
    parser.add_argument(
        '--run_id',
        help = 'Long ID for the MiSeq run',
        metavar = '<run_id>',
        dest = 'run_id',
        required = True
    )
    
    parser.add_argument(
        '--expt_id',
        help = 'Experiment ID for the run',
        metavar = '<expt_id>',
        dest = 'expt_id',
        required = True
    )
    
    parser.add_argument(
        '--panel',
        help = 'Type of panel',
        metavar = '<panel>',
        dest = 'panel',
        choices = ['myeloid', 'lymphoid_CLL_MZL'],
        required = True
    )
        
    args = parser.parse_args()
    
    return args


def call_subprocess(cmd):
    
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    with p.stdout:
        for line in iter(p.stdout.readline, b''):
            logger.info(line.replace('\n', ''))
    
    p.wait()
    
    return None


class Pipeline(object):
    
    def __init__(self, args):
        
        # core attributes
        self.run_id = args.run_id
        self.expt_id = args.expt_id
        self.panel = args.panel
                
        # set the run and experiment directories
        if self.panel == 'myeloid':
            
            self.raw_data_dir = os.path.join(BACKUP_DIR, 'Myeloid_output', self.run_id)
            self.run_dir = os.path.join(ANALYSIS_DIR, 'myeloid', self.expt_id)
            
        elif self.panel == 'lymphoid_CLL_MZL':
            
            self.raw_data_dir = os.path.join(BACKUP_DIR, 'Lymphoid_output', self.run_id)
            self.run_dir = os.path.join(ANALYSIS_DIR, 'lymphoid', self.expt_id)

        # set the reference files directory (variable only needed temporarily)
        refs = os.path.join(REFERENCES_DIR, self.panel)
        
        # temporary function to grab some file paths
        def fn(x):
            return os.path.join(
                refs, next(f for f in os.listdir(refs) if self.panel + x in f)
            )
        
        # set file paths
        self.targets_file = fn('_targets_1_H')
        self.manifest_file = fn('_manifest_0_nH')
        self.gene_list = fn('_gene_list')
        self.trans_list = fn('_transcripts')
        self.hotspots_file = fn('_hotspots')
        self.exclusions_file = fn('_exclusions')
        
        # sub-directory paths
        self.alamut_dir = os.path.join(self.run_dir, 'alamut')
        self.bam_dir = os.path.join(self.run_dir, 'bam')
        self.coverage_dir = os.path.join(self.run_dir, 'coverage')
        self.qc_dir = os.path.join(self.run_dir, 'QC')
        self.results_dir = os.path.join(self.run_dir, 'results')
        self.vcf_dir = os.path.join(self.run_dir, 'vcf')
        
        self.subdir_list = [
            self.alamut_dir,
            self.bam_dir,
            self.coverage_dir,
            self.qc_dir,
            self.results_dir,
            self.vcf_dir
        ]
        
    
    def run(self):
        
        start = strftime("%Y-%m-%d %H:%M:%S")
        user = getpass.getuser()
        
        logger.info('Run started at %s by %s' % (start, user))
        
        self.check_files()
        self.setup()
        self.run_fastqc()
        self.run_bwa()
        self.run_samtools()
        self.run_picard()
        self.run_varscan()
        
        if self.panel == 'myeloid':
            self.asxl1_check()
        
        self.run_filtervcf()
        self.run_alamut()
        self.generate_results()
        
        end = strftime("%Y-%m-%d %H:%M:%S")
        
        logger.info('Run started at %s' % start)
        logger.info('Run ended at %s' % end)
        
        log_fname = self.expt_id + '.log'
        
        shutil.copyfile(
            os.path.join(ANALYSIS_DIR, 'pipeline_log', log_fname),
            os.path.join(self.run_dir, log_fname)
        )
        
    
    def check_files(self):
        
        '''Check that the required files and directories exist'''
        
        # check that the raw data directory exists
        if os.path.isdir(self.raw_data_dir):
            logger.info('Raw data directory: %s' % self.raw_data_dir)
        else:
            logger.error('Raw data directory does not exist: %s' % self.raw_data_dir)
            sys.exit()

        # check that the run ID folder does not already exist
        if not os.path.isdir(self.run_dir):
            os.mkdir(self.run_dir)
            logger.info('Created run directory: %s' % self.run_dir)
        else:
            logger.error('Run directory already exists: %s' % self.run_dir)
            sys.exit()

        # check that the targets and manifest files exist
        files = {
            'Targets file': self.targets_file,
            'Manifest file': self.manifest_file,
            'Gene list file': self.gene_list,
            'Transcripts file': self.trans_list,
            'Hotspots file': self.hotspots_file,
            'Exclusions file': self.exclusions_file
        }

        for k, v in files.iteritems():
            if os.path.exists(v):
                logger.info('%s: %s' % (k, v))
            else:
                logger.error('%s does not exist: %s' % (k, v))
                sys.exit()
        
        return None
    
    
    def setup(self):
        
        '''Set up the directory structure and copy SampleSheet.csv'''
        
        # cd to the run directory
        os.chdir(self.run_dir)

        # set up the directory structure
        logger.info('Setting up the directory structure')
        for subdir in self.subdir_list: os.mkdir(subdir)
        
        # copy SampleSheet.csv to the run directory
        logger.info('Copying SampleSheet.csv to %s' % self.run_dir)
        shutil.copyfile(
            os.path.join(self.raw_data_dir, 'SampleSheet.csv'),
            os.path.join(self.run_dir, 'SampleSheet.csv')
        )

        # get a list of samples from SampleSheet
        sample_sheet = os.path.join(self.run_dir, 'SampleSheet.csv')
        
        # work out what line the string '[Data]' appears on in SampleSheet
        n = next(i for i, line in enumerate(open(sample_sheet, 'r')) if '[Data]' in line)
        
        # skip n lines before loading SampleSheet data into a df
        df = pd.read_csv(sample_sheet, skiprows=n+1, usecols=['Sample_ID'])
        
        # replace '/' with '-' in Sample_ID column and send to list
        self.sample_list = df['Sample_ID'].str.replace('/', '-').tolist()

        logger.info('Samples to be processed:')
        
        for s in self.sample_list:
            logger.info(s)

        return None
    
    
    def run_fastqc(self):
        
        '''Runs FastQC'''

        logger.info('>>>FastQC>>>')

        # make the fastqc_reports directory
        fastqc_reports_dir = os.path.join(self.qc_dir, 'fastqc_reports')
        os.mkdir(fastqc_reports_dir)

        # cd to the fastq data directory
        fastq_dir = os.path.join(self.raw_data_dir, 'Data', 'Intensities', 'BaseCalls')
        os.chdir(fastq_dir)

        # get a list of fastq files to iterate over
        self.fastq_files = [f for f in os.listdir(fastq_dir) if '.fastq.gz' in f]

        l = len(self.fastq_files)
        
        for c, f in enumerate(self.fastq_files):
            
            logger.info('Running FastQC for file: %s (%d of %d)' % (f, c+1, l))
            
            # -t 8 specifies the number of threads
            cmd = [FASTQC, f, '-o', fastqc_reports_dir, '-t', '8']
            
            call_subprocess(cmd)

        return None
    
    
    def run_bwa(self):
        
        '''Runs BWA MEM'''
        
        logger.info('>>>BWA MEM>>>')

        R1_list = [i for i in self.fastq_files if '_R1_' in i]
        R2_list = [i for i in self.fastq_files if '_R2_' in i]

        self.sam_file_list = []

        l = len(self.sample_list)
        
        for c, sample in enumerate(self.sample_list):

            logger.info('Aligning sample %s (%d of %d)' % (sample, c+1, l))

            fastq_R1 = next(i for i in R1_list if i.startswith(sample))
            fastq_R2 = next(i for i in R2_list if i.startswith(sample))

            sam_file = sample + '.sam'
            
            logger.info('Aligning files: %s and %s' % (fastq_R1, fastq_R2))
            logger.info('Reference genome: %s' % REF_GENOME)
            
            cmd = ' '.join(
                [
                    BWA, 'mem',
                    '-t', '8',
                    REF_GENOME,
                    fastq_R1, fastq_R2,
                    '>', os.path.join(self.bam_dir, sam_file)
                ]
            )
            
            #####################################################################
            #
            # Note that Shell = True for the following command. This is due to
            # the use of '>' redirection at the command line.
            #
            # Note also that when using Shell=True, subprocess.Popen expects
            # a string and not a list.
            #
            # from the docs:
            #
            # "The shell argument (which defaults to False) specifies whether
            # to use the shell as the program to execute. If shell is True, it
            # is recommended to pass args as a string rather than as a sequence."
            #
            # https://docs.python.org/2/library/subprocess.html#popen-constructor
            #
            #####################################################################
            
            proc = subprocess.Popen(
                cmd,
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
                shell = True
            )
                    
            proc_output, _ = proc.communicate()
            logger.info(proc_output)
            logger.info('Output SAM file: %s' % sam_file)

            self.sam_file_list.append(sam_file)

        return None
    
    
    def run_samtools(self):
        
        '''
        SAMtools view: SAM to BAM conversion
        SAMtools sort: creates sorted BAM file
        SAMtools index: creates BAM index file
        SAMtools mpileup: creates mpileup file
        '''
        
        logger.info('>>>SAMtools>>>')
        
        os.chdir(self.bam_dir)
    
        self.bam_list = []
        self.pileup_list = []
        
        l = len(self.sam_file_list)

        for c, f in enumerate(self.sam_file_list):

            logger.info('Processing %s (%d of %d)' % (f, c+1, l))

            # set the file names
            bam_file = f.replace('.sam', '.bam')
            sorted_bam_file = f.replace('.sam', '_sorted.bam')
            pileup_file = f.replace('.sam', '.pileup')
            pileup_path = os.path.join(self.coverage_dir, pileup_file)

            # run samtools view
            logger.info('Converting %s to %s' % (f, bam_file))
            cmd = [SAMTOOLS, 'view', '-b', '-@', '8', f, '-o', bam_file]
            call_subprocess(cmd)

            # run samtools sort
            logger.info('Sorting %s' % bam_file)
            cmd = [SAMTOOLS, 'sort', '-@', '8', bam_file, '-o', sorted_bam_file]
            call_subprocess(cmd)

            # run samtools index
            logger.info('Indexing %s' % sorted_bam_file)
            cmd = [SAMTOOLS, 'index', sorted_bam_file]
            call_subprocess(cmd)
            
            # run samtools mpileup on the sorted bam file
            logger.info('Creating pileup file from %s' % sorted_bam_file)
            cmd = [
                SAMTOOLS, 'mpileup', '-B',
                '-f', REF_GENOME,
                '-l', self.manifest_file,
                sorted_bam_file,
                '-o', pileup_path
            ]
            call_subprocess(cmd)
            
            self.bam_list.append(sorted_bam_file)
            self.pileup_list.append(pileup_file)

        return None
    
    
    def run_picard(self):
        
        '''Runs Picard CollectAlignmentSummaryMetrics and CollectHsMetrics'''
        
        logger.info('>>>Picard>>>')

        mapping_reports_dir = os.path.join(self.bam_dir, 'mapping_reports')
        os.mkdir(mapping_reports_dir)

        l = len(self.bam_list)
        
        for c, f in enumerate(self.bam_list):

            # run CollectAlignmentSummaryMetrics
            logger.info('Running CollectAlignmentSummaryMetrics for %s (%d of %d)' % (f, c+1, l))

            align_sum_file = f.replace('_sorted.bam', '_alignSummary.txt')
            align_sum_path = os.path.join(mapping_reports_dir, align_sum_file)

            cmd = [
                'java', '-Xmx8g', '-jar', PICARD,
                'CollectAlignmentSummaryMetrics',
                'I=', f,
                'O=', align_sum_path
            ]

            call_subprocess(cmd)
            
            logger.info('AlignmentSummaryMetrics output saved to %s' % align_sum_file)

            # run CollectHsMetrics
            logger.info('Running CollectHsMetrics for %s' % f)

            targets_sum_file = f.replace('_sorted.bam', '_targetsSummary.txt')
            targets_sum_path = os.path.join(mapping_reports_dir, targets_sum_file)

            cmd = [
                'java', '-Xmx8g', '-jar', PICARD,
                'CollectHsMetrics',
                'I=', f,
                'O=', targets_sum_path,
                'R=', REF_GENOME,
                'BAIT_INTERVALS=', self.targets_file,
                'TARGET_INTERVALS=', self.targets_file
            ]
            
            call_subprocess(cmd)
            
            logger.info('HsMetrics output saved to %s' % targets_sum_file)

        return None
    
    
    def run_varscan(self):
        
        '''Runs VarScan2'''
        
        logger.info('>>>VarScan2>>>')

        os.chdir(self.coverage_dir)

        self.vcf_list = []

        l = len(self.pileup_list)
        
        for c, f in enumerate(self.pileup_list):
            
            logger.info('Running VarScan2 for %s (%d of %d)' % (f, c+1, l))

            vcf_file = f.replace('.pileup', '.vcf')
            vcf_file_path = os.path.join(self.vcf_dir, vcf_file)

            cmd = ' '.join(
                [
                    'java', '-Xmx8g', '-jar', VARSCAN,
                    'mpileup2cns', f,
                    '--min-coverage', '50',
                    '--min-reads2', '5',
                    '--min-var-freq', '0.02',
                    '--p-value', '0.98',
                    '--strand-filter', '0',
                    '--variants', '1',
                    '--output-vcf', '1',
                    '>', vcf_file_path
                ]
            )
            
            # as with BWA, Varscan has to be run with Shell=True
            # (so don't use the call_subprocess function)
            proc = subprocess.Popen(
                cmd,
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
                shell = True
            )
                    
            proc_output, _ = proc.communicate()
            logger.info(proc_output)
            logger.info('VCF file saved to %s' % vcf_file)
            
            self.vcf_list.append(vcf_file)

        return None
    
    
    def asxl1_check(self):
        
        '''Checks for the ASXL1 23bp deletion'''
        
        logger.info('>>>Checking for ASXL1 23bp deletion>>>')
        
        os.chdir(self.vcf_dir)
        
        realignment_list = []
        
        for vcf in self.vcf_list:
            
            logger.info('Checking %s' % vcf)
            
            with open(vcf, 'r') as f:
                lines = [l for l in f.readlines() if "##" in l]
                
            df = pd.read_table(vcf, sep='\t', skiprows=len(lines))
            
            df_asxl1_missense = df.loc[
                (df['POS'] == 31022415) &
                (df['REF'] == 'A') &
                (df['ALT'] == 'C')
            ]

            if df_asxl1_missense.empty:
                logger.info('ASXL1 31022415_A/C not found in %s' % vcf)
            else:
                logger.info('ASXL1 31022415_A/C found in %s' % vcf)
                realignment_list.append(vcf.replace('.vcf', ''))
                
        if realignment_list != []:
            
            realigned_dir = os.path.join(self.run_dir, 'realigned')
            
            if not os.path.exists(realigned_dir): os.mkdir(realigned_dir)
                
            for sample in realignment_list:

                logger.info('Realigning %s' % sample)

                bam_abra_in = sample + '_sorted.bam'
                bam_abra_out = sample + '_realigned.bam'
                
                # run ABRA
                logger.info('Running ABRA on %s' % sample)
                
                cmd = [
                    'java', '-Xmx8g', '-jar', ABRA,
                    '--in', os.path.join(self.run_dir, 'bam', bam_abra_in),
                    '--out', os.path.join(self.run_dir, 'realigned', bam_abra_out),
                    '--ref', REF_GENOME,
                    '--dist', '150',
                    '--targets', ASXL1_BED,
                    '--threads', '8'
                ]
                call_subprocess(cmd)
                
                os.chdir(realigned_dir)
                
                # run samtools sort
                logger.info('Sorting %s' % bam_abra_out)
                
                realigned_sorted_bam = bam_abra_out.replace('.bam', '_sorted.bam')
                
                cmd = [
                    SAMTOOLS, 'sort',
                    '-@', '8',
                    bam_abra_out,
                    '-o', realigned_sorted_bam
                ]
                call_subprocess(cmd)

                # run samtools index
                logger.info('Indexing %s' % realigned_sorted_bam)
                
                cmd = [SAMTOOLS, 'index', realigned_sorted_bam]
                
                call_subprocess(cmd)

                # run samtools mpileup on the sorted bam file
                logger.info('Creating pileup file from %s' % realigned_sorted_bam)
                
                realigned_pileup = sample + '_realigned.pileup'
                
                cmd = [
                    SAMTOOLS, 'mpileup', '-B',
                    '-f', REF_GENOME,
                    '-l', self.manifest_file,
                    realigned_sorted_bam,
                    '-o', realigned_pileup
                ]
                call_subprocess(cmd)
                
                logger.info('Running VarScan2 on %s' % realigned_pileup)
                
                realigned_vcf = sample + '_realigned.vcf'
                
                cmd = ' '.join(
                    [
                        'java', '-Xmx8g', '-jar', VARSCAN,
                        'mpileup2cns', realigned_pileup,
                        '--min-coverage', '50',
                        '--min-reads2', '5',
                        '--min-var-freq', '0.02',
                        '--p-value', '0.98',
                        '--strand-filter', '0',
                        '--variants', '1',
                        '--output-vcf', '1',
                        '>', realigned_vcf
                    ]
                )

                # VarScan must be run with Shell=True, so don't use call_subprocess
                proc = subprocess.Popen(
                    cmd,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.STDOUT,
                    shell = True
                )

                proc_output, _ =  proc.communicate()
                logger.info(proc_output)
                
                # look for the 23bp deletion in the realigned vcf
                with open(realigned_vcf, 'r') as f:
                    lines = [l for l in f.readlines() if "##" in l]

                df = pd.read_table(realigned_vcf, sep='\t', skiprows=len(lines))
                
                df_asxl1_23bp_del = df.loc[
                    (df['POS'] == 31022402) &
                    (df['REF'] == 'TCACCACTGCCATAGAGAGGCGGC') &
                    (df['ALT'] == 'T')
                ]

                if df_asxl1_23bp_del.empty:
                    logger.info('ASXL1 23bp deletion not found in %s' % realigned_vcf)
                    
                else:
                    logger.info('ASXL1 23bp deletion found in %s' % realigned_vcf)
                    
                    original_vcf = os.path.join(self.run_dir, 'vcf', sample + '.vcf')
                    new_vcf = os.path.join(self.run_dir, 'realigned', sample + '.vcf')
                    
                    # read in the original vcf file but only keep the header lines
                    with open(original_vcf, 'r') as f:
                        lines = [l for l in f.readlines() if "##" in l]

                    # write the header lines to a new file
                    with open(new_vcf, 'w') as f:
                        f.writelines(lines)

                    # read only the VCF data rows into a df by skipping the header rows
                    df_orig = pd.read_table(original_vcf, sep='\t', skiprows=len(lines))
                    
                    # drop the asxl1 missense variant - NOTE THE ~TILDE~ MEANING NEGATION
                    df_orig = df_orig.loc[
                        ~(
                            (df_orig['POS'] == 31022415) &
                            (df_orig['REF'] == 'A') &
                            (df_orig['ALT'] == 'C')
                        )
                    ]
                    
                    # append the row relating to the 23bp deletion to the original vcf
                    df_out = df_orig.append(df_asxl1_23bp_del, ignore_index=True)
                    
                    # adjust format of CHROM for sorting
                    df_out['#CHROM'] = df_out['#CHROM'].str.replace('chr', '')
                    df_out['#CHROM'] = df_out['#CHROM'].str.replace('X', '23')
                    df_out['#CHROM'] = pd.to_numeric(df_out['#CHROM'])

                    # sort on CHROM and POS
                    df_out.sort_values(['#CHROM', 'POS'], axis=0, inplace=True)

                    # revert to original format of CHROM before writing to file
                    df_out['#CHROM'] = df_out['#CHROM'].apply(str)
                    df_out['#CHROM'] = df_out['#CHROM'].str.replace('23', 'X')
                    df_out['#CHROM'] = 'chr' + df_out['#CHROM']

                    # write to file
                    df_out.to_csv(new_vcf, sep='\t', index=False, mode='a')
                    
                    # copy the original vcf to realigned dir
                    shutil.copyfile(
                        original_vcf,
                        os.path.join(realigned_dir, sample + '_original.vcf')
                    )
                    
                    # copy the new vcf to vcf dir (overwriting the original)
                    shutil.copyfile(new_vcf, original_vcf)
                    
        return None
    
    
    def run_filtervcf(self):
        
        '''Runs filtervcf.py'''
        
        logger.info('>>>VCF Filtering>>>')

        os.chdir(self.vcf_dir)
        os.mkdir('filtered')

        for f in self.vcf_list:
            
            logger.info('Filtering %s' % f)
            
            vcf_out = os.path.join(self.vcf_dir, 'filtered', f)
            
            cmd = [
                'python', FILTERVCF,
                '--in', f,
                '--out', vcf_out,
                '--hotspots', self.hotspots_file,
                '--exclusions', self.exclusions_file
            ]
            
            call_subprocess(cmd)

        return None
        
    
    def run_alamut(self):
        
        '''Runs Alamut Batch'''
        
        logger.info('>>>Alamut Batch>>>')

        os.chdir(os.path.join(self.vcf_dir, 'filtered'))
        
        ann_dir = os.path.join(self.alamut_dir, 'ann')
        os.mkdir(ann_dir)
        
        unann_dir = os.path.join(self.alamut_dir, 'unann')
        os.mkdir(unann_dir)

        l = len(self.vcf_list)
        
        for c, f in enumerate(self.vcf_list):
            
            logger.info('Running Alamut Batch for %s (%d of %d)' % (f, c+1, l))

            ann_file = f.replace('.vcf', '.ann')
            ann_path = os.path.join(ann_dir, ann_file)
            
            unann_file = f.replace('.vcf', '.unann')
            unann_path = os.path.join(unann_dir, unann_file)
            
            cmd = [
                ALAMUT,
                '--in', f,
                '--ann', ann_path,
                '--unann', unann_path,
                '--glist', self.gene_list,
                '--translist', self.trans_list,
                '--processes', '8'
            ]
            
            call_subprocess(cmd)

        return None
    
    
    def generate_results(self):
        
        '''Generate results'''
        
        logging.info('>>>Results>>>')

        for f in self.vcf_list:

            logger.info('Generating results for sample %s' % f.replace('.vcf', ''))

            ann_file = f.replace('.vcf', '.ann')
            ann_path = os.path.join(self.alamut_dir, 'ann', ann_file)
            
            out_file = f.replace('.vcf', '.txt')
            out_path = os.path.join(self.results_dir, out_file)

            cmd = [
                'python', GENERATE_RESULTS,
                '--vcf', f,
                '--ann', ann_path,
                '--out', out_path,
                '--hotspots', self.hotspots_file
            ]
                        
            call_subprocess(cmd)

        return None
    
    
if __name__ == "__main__":
    args = get_args()
    
    # start logging
    logFormatter = logging.Formatter(
        fmt='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    logfile = os.path.join(ANALYSIS_DIR, 'pipeline_log', args.expt_id + '.log')
    
    fileHandler = logging.FileHandler(logfile, mode='w')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    
    # run the pipeline
    pipeline = Pipeline(args)
    pipeline.run()