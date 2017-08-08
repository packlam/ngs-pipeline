import argparse
import os
import sys
import pandas as pd


# directories
REFERENCES_DIR = '/HMDS_BIOINFORMATICS/references/'


def get_args():
    
    parser = argparse.ArgumentParser(
        description = 'Filter a VCF file'
    )
    
    parser.add_argument(
        '--in',
        help = 'Path to input VCF file',
        metavar = '<input.vcf>',
        dest = 'vcf_in',
        required = True
    )
    
    parser.add_argument(
        '--out',
        help = 'Path to output VCF file',
        metavar = '<output.vcf>',
        dest = 'vcf_out',
        required = True
    )
    
    parser.add_argument(
        '--hotspots',
        help = 'Hotspots file',
        metavar = '<hotspots>',
        dest = 'hotspots',
        required = True
    )
    
    parser.add_argument(
        '--exclusions',
        help = 'Exclusions file',
        metavar = '<exclusions>',
        dest = 'exclusions',
        required = True
    )
    
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError('VAF must be between 0.0 and 1.0')
        return x
    
    parser.add_argument(
        '--vaf',
        help = 'Minimum allele frequency (default: 0.05)',
        metavar = '<vaf>',
        dest = 'vaf',
        type = restricted_float,
        default = 0.05,
        required = False
    )
    
    parser.add_argument(
        '--depth',
        help = 'Minimum read depth (default: 100)',
        metavar = '<depth>',
        dest = 'depth',
        type = int,
        default = 100,
        required = False
    )
    
    args = parser.parse_args()
    
    return args


def load_from_txt(file):
    
    '''Load hotspot or exclusion files into a list'''
    
    with open(file, 'r') as f:
        data = f.read().splitlines()
        data = [int(i) for i in data] # convert str to int
    
    return data


class Filter(object):
    
    def __init__(self, args):
        
        self.vcf_in = args.vcf_in
        self.vcf_out = args.vcf_out
        self.vcf_out_dir = os.path.dirname(self.vcf_out)
        self.hotspots_file = args.hotspots
        self.exclusions_file = args.exclusions
        
        # filter values
        self.vaf_filter = args.vaf
        self.depth_filter = args.depth
        
        return None
    
    
    def run(self):
        
        # run checks on required directories and files
        self.check_files()
        
        # load data
        df_split = self.split_vcf()
        df_filt = self.filter_variants(df_split)
        df_hotspots = self.extract_hotspots(df_split)
        df_out = self.concatenate_dfs(df_filt, df_hotspots)
        
        # write in tab separated format to vcf_out in append mode
        print('\nSaving VCF file: %s' % self.vcf_out)
        df_out.to_csv(self.vcf_out, sep='\t', index=False, mode='a')
        
        return None
    
    
    def check_files(self):
        
        files = {
            'Input VCF file': self.vcf_in,
            'Hotspots file': self.hotspots_file,
            'Exclusions file': self.exclusions_file
        }

        # check that input files exist
        for k, v in files.iteritems():
            if os.path.isfile(v):
                print('\nOK: %s exists: %s' % (k, v))
            else:
                sys.exit('\nERROR: %s does not exist: %s' % (k, v))

        # check that the specified output directory exists
        if os.path.isdir(self.vcf_out_dir):
            print('\nOK: Output VCF directory exists: %s' % self.vcf_out_dir)
        else:
            sys.exit('\nERROR: Output VCF directory does not exist: %s' % self.vcf_out_dir)

        # check that the specified input and output have .vcf suffix
        for f in [self.vcf_in, self.vcf_out]:
            if not f.endswith('.vcf'):
                sys.exit('File must have a .vcf suffix: %s' % f)
        
        return None
    
    
    def split_vcf(self):
        
        '''Splits the data in the 'Sample1' column into separate columns'''
        
        df_vcf = self.load_vcf()
        
        # generate list of column names by splitting 'FORMAT' column on ':' delimiter
        self.format_labels = df_vcf.loc[0, 'FORMAT'].split(':')

        # create new df by splitting out the 'Sample1' column on ':' delimiter
        df_split = df_vcf['Sample1'].str.split(':', expand=True)

        # apply column names to df_split
        df_split.columns = self.format_labels

        # convert FREQ % values to decimal (e.g. 25% to 0.25)
        df_split['FREQ'] = df_split['FREQ'].str.replace('%', '') # remove % symbols
        df_split['FREQ'] = pd.to_numeric(df_split['FREQ']) # convert str to float
        df_split['FREQ'] = df_split['FREQ'] / 100 # divide by 100

        # convert DP values to int
        df_split['DP'] = pd.to_numeric(df_split['DP']) # convert str to int

        # drop the original 'Sample1' column from df_vcf
        df_vcf = df_vcf.drop('Sample1', 1)

        # concatenate df_vcf and df_split side by side
        df = pd.concat([df_vcf, df_split], axis=1)

        return df
    
    
    def load_vcf(self):
        
        '''Write the VCF header to vcf_out and read the data lines into a df'''

        # read in the VCF file but only keep the header lines
        with open(self.vcf_in, 'r') as f:
            lines = [l for l in f.readlines() if "##" in l]

        # write the header lines to a new file
        with open(self.vcf_out, 'w') as f:
            f.writelines(lines)

        # read only the VCF data rows into a df by skipping the header rows
        df = pd.read_table(self.vcf_in, sep='\t', skiprows=len(lines))
        
        return df
    
    
    def filter_variants(self, df_split):
        
        '''Filters the variants on VAF and depth, and removes excluded variants'''
        
        # load exclusions file
        exclusions = load_from_txt(self.exclusions_file)

        # remove any variants in the exclusions list
        df = df_split[~df_split['POS'].isin(exclusions)]

        # keep only variants with VAF >= 0.05
        df = df[df['FREQ'] >= self.vaf_filter]

        # keep only variants with reads >= 100
        df = df[df['DP'] >= self.depth_filter]

        return df
    
    
    def extract_hotspots(self, df_split):
        
        '''Returns a df containing only variants in the hotspot list'''

        # load the hotspots file
        hotspots = load_from_txt(self.hotspots_file)

        # df = all the rows in df_split where the position is in hotspots    
        df = df_split[df_split['POS'].isin(hotspots)]

        return df
    
    
    def concatenate_dfs(self, df_filt, df_hotspots):
        
        # add df_hotspots to the bottom of df_filt
        df = pd.concat([df_filt, df_hotspots], axis=0, ignore_index=True)

        # drop any duplicate rows
        df.drop_duplicates(inplace=True)

        # convert 'FREQ' from float to str in order to carry out next step
        df['FREQ'] = df['FREQ'].astype(str)

        # convert 'DP' from int to str in order to carry out next step
        df['DP'] = df['DP'].astype(str)

        # recreate the 'Sample1' column by sticking the constituent columns together
        df['Sample1'] = df[self.format_labels].apply(lambda x: ':'.join(x), axis=1)

        # drop the unwanted columns
        df = df.drop(self.format_labels, 1)
        
        # set the dtype of the POS column to string
        df['POS'] = df['POS'].astype('str')

        # specify the columns with which to create the ID
        cols = ['#CHROM', 'POS', 'REF', 'ALT']

        # create an ID for each variant by concatenating chrom, pos, ref and alt
        df['ID'] = df[cols].apply(lambda x: '_'.join(x), axis=1)
        
        return df
    
    
if __name__ == "__main__":
    args = get_args()
    filt = Filter(args)
    filt.run()