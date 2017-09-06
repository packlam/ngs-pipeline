import argparse
import os
import sys
from collections import Counter
import numpy as np
import pandas as pd


def get_args():
    
    parser = argparse.ArgumentParser(
        description = 'Generate results from a VCF file and Alamut Batch results file'
    )
    
    parser.add_argument(
        '--vcf',
        help = 'Path to input VCF file',
        metavar = '<input.vcf>',
        dest = 'vcf_in',
        required = True
    )
    
    parser.add_argument(
        '--ann',
        help = 'Path to input ann file',
        metavar = '<input.ann>',
        dest = 'ann_in',
        required = True
    )
    
    parser.add_argument(
        '--out',
        help = 'Path to output file',
        metavar = '<output.txt>',
        dest = 'results_out',
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
    
    args = parser.parse_args()
    
    return args


class ResultsGenerator(object):
    
    def __init__(self, args):
        
        self.vcf_in = args.vcf_in
        self.ann_in = args.ann_in
        self.results_out = args.results_out
        self.results_dir = os.path.dirname(self.results_out)
        self.hotspots_file = args.hotspots
        self.exclusions_file = args.exclusions
                
        return None
    
    
    def run(self):
        
        # check that require files and directories exist
        self.check_files()
        
        # load the vcf and alamut batch data
        df_vcf = self.load_vcf()
        df_ann = self.load_ann()
        
        # concatenate df_vcf and df_ann side by side
        df_merged = pd.merge(
            df_vcf, df_ann,
            how = 'inner',
            on = ['id', 'chr', 'POS', 'REF', 'ALT']
        )
        
        # create a new column named 'keep', assign all values to False
        df_merged['keep'] = False
        
        # discard variants on the exclusions list
        df_merged = self.discard_exclusions(df_merged)
        
        # retain variants on the hotspots list
        df_merged = self.retain_hotspots(df_merged)

        # run the df through all the functions
        df_merged = self.cosmic(df_merged)
        df_merged = self.splice(df_merged)
        df_merged = self.filter_variants(df_merged)      
        
        if df_merged.empty:
            
            print('\nNo variants remain after filtering')
            
        else:
        
            # tweak the output
            df_merged = self.tweak(df_merged)

            # save the output
            print('\nSaving results file: %s' % self.results_out)
            df_merged.to_csv(self.results_out, sep='\t', header=True, index=False)
        
        return None
    
    
    def check_files(self):
        
        '''Checks that the required files and directories exist'''
        
        files = {
            'Input VCF file': self.vcf_in,
            'Input Alamut Batch file': self.ann_in,
            'Hotspots file': self.hotspots_file,
            'Exclusions file': self.exclusions_file,
        }

        # check that input files exist
        for k, v in files.iteritems():
            if os.path.isfile(v):
                print('\nOK: %s exists: %s' % (k, v))
            else:
                sys.exit('\nERROR: %s does not exist: %s' % (k, v))

        # check that the specified output directory exists
        if os.path.isdir(self.results_dir):
            print('\nOK: Results output directory exists: %s' % self.results_dir)
        else:
            sys.exit('\nERROR: Results output directory does not exist: %s' % self.results_dir)

        # check that all file paths end with the correct extension
        extensions = {self.vcf_in: '.vcf', self.ann_in: '.ann', self.results_out: '.txt'}
        for file, ext in extensions.iteritems():
            if not file.endswith(ext):
                sys.exit('File must have a %s suffix: %s' % (ext, file))
        
        return None
        
    
    def load_vcf(self):
        
        '''Loads the VCF data into a dataframe'''
        
        print('\nLoading VCF data')

        # read in the VCF file but only keep the header lines
        with open(self.vcf_in, 'r') as f:
            lines = [l for l in f.readlines() if "##" in l]
            numheadlines = len(lines)

        df_vcf = pd.read_table(self.vcf_in, sep='\t', skiprows=numheadlines)

        # generate list of column names by splitting 'FORMAT' column on ':' delimiter
        format_labels = df_vcf.loc[0, 'FORMAT'].split(':')

        # create new df by splitting out the 'Sample1' column on ':' delimiter
        df_split = df_vcf['Sample1'].copy().str.split(':', expand=True)

        # apply column names to df_split
        df_split.columns = format_labels
        
        # convert FREQ % values to decimal (e.g. 25% to 0.25)
        df_split['FREQ'] = df_split['FREQ'].str.replace('%', '') # remove % symbols
        df_split['FREQ'] = pd.to_numeric(df_split['FREQ']) # convert str to float
        df_split['FREQ'] = df_split['FREQ'] / 100 # divide by 100
        df_split['FREQ'] = df_split['FREQ'].round(decimals=3)

        # convert 'DP' from str to int
        df_split['DP'] = pd.to_numeric(df_split['DP'])

        # concatenate df_vcf and df_split side by side
        df = pd.concat([df_vcf, df_split], axis=1)

        # assign the sample id
        sample_id = os.path.basename(self.vcf_in)
        
        # modify the sample_id as necessary
        for x, y in {'.vcf': '', '-': '/', '_f': ''}.iteritems():
            sample_id = sample_id.replace(x, y)
        
        # add sample id column
        df['sample_id'] = sample_id
        
        # set the dtype of the POS column to string
        df['POS'] = df['POS'].astype('str')

        # specify the columns with which to create the ID
        cols = ['#CHROM', 'POS', 'REF', 'ALT']

        # create an ID for each variant by concatenating chrom, pos, ref and alt
        df['ID'] = df[cols].apply(lambda x: '_'.join(x), axis=1)
        
        # set the dtype of the POS column back to int
        df['POS'] = df['POS'].astype('int')
        
        # rename columns
        df.rename(columns={'#CHROM': 'chr', 'ID': 'id'}, inplace=True)

        # keep only desired columns
        df = df[['sample_id', 'id', 'chr', 'POS', 'REF', 'ALT', 'DP', 'FREQ']]

        return df
    
    
    def load_ann(self):
        
        '''Loads a subset of the Alamut Batch data into a dataframe'''
        
        print('\nLoading Alamut Batch data')

        cols = [
            'chrom',
            'inputPos',
            'inputRef',
            'inputAlt',
            'gene',
            'varType',
            'codingEffect',
            'cNomen',
            'pNomen',
            'rsId',
            'cosmicIds',
            'cosmicTissues',
            'cosmicFreqs',
            'cosmicSampleCounts',
            'distNearestSS',
            'nearestSSChange',
            '1000g_AF',
            'exacAltFreq_all'
        ]

        # load the alamut data into a df
        df = pd.read_table(self.ann_in, sep='\t', usecols=cols)

        # rename columns
        df.rename(
            columns={
                'chrom': 'chr',
                'inputPos': 'POS',
                'inputRef': 'REF',
                'inputAlt': 'ALT',
                'codingEffect': 'consequence'
            },
            inplace=True
        )

        # remove spaces from some of the possible values in 'consequence'
        df.replace(
            {
                'stop gain': 'stop_gain',
                'start loss': 'start_loss',
                'stop loss': 'stop_loss'
            },
            inplace=True
        )
        
        # convert chr values from e.g. 1 --> chr1
        df['chr'] = df['chr'].astype(str)
        df['chr'] = 'chr' + df['chr']
        
        # set the dtype of the POS column to string
        df['POS'] = df['POS'].astype('str')

        # specify the columns with which to create the ID
        cols = ['chr', 'POS', 'REF', 'ALT']

        # create an ID for each variant by concatenating chrom, pos, ref and alt
        df['id'] = df[cols].apply(lambda x: '_'.join(x), axis=1)
        
        # set the dtype of the POS column back to int
        df['POS'] = df['POS'].astype('int')

        return df
    
    
    def discard_exclusions(self, df_merged):
        
        '''Removes variants that are in the exclusions file'''
        
        cols = ['chr', 'POS', 'REF', 'ALT']
        
        df_exc = pd.read_table(self.exclusions_file, sep='\t', usecols=cols)
        
        df = pd.merge(
            df_merged, df_exc,
            on = ['chr', 'POS', 'REF', 'ALT'],
            how = 'left',
            indicator = 'excluded'
        )
        
        df['excluded'] = np.where(df['excluded'] == 'both', True, False)
        
        df = df.loc[df['excluded'] == False]
        
        df.drop('excluded', axis=1, inplace=True)
        
        return df
    
    
    def retain_hotspots(self, df_merged):
        
        '''Retains the variants that are in the hotspots file'''
        
        cols = ['chr', 'POS', 'REF', 'ALT']
        
        df_hot = pd.read_table(self.hotspots_file, sep='\t', usecols=cols)
        
        df = pd.merge(
            df_merged, df_hot,
            on = ['chr', 'POS', 'REF', 'ALT'],
            how = 'left',
            indicator = 'hotspot'
        )
        
        df['hotspot'] = np.where(df['hotspot'] == 'both', True, False)
        
        df.loc[
            (df['hotspot'] == True) &
            (df['FREQ'] >= 0.02),
            'keep'
        ] = True
        
        df.drop('hotspot', axis=1, inplace=True)
        
        return df
        
    
    def cosmic(self, df_whole):
        
        '''Does some clever magic with the COSMIC data'''
        
        # specify the columns you want to keep
        cols = ['id', 'cosmicIds', 'cosmicTissues',
                'cosmicFreqs', 'cosmicSampleCounts']

        # keep only these columns
        df = df_whole[cols].copy()

        # keep only the rows that have at least one cosmic id
        df = df[pd.notnull(df['cosmicIds'])]

        # convert all values in these columns to lists by splitting on the '|' character        
        df['cosmicIds'] = df['cosmicIds'].str.split('|').tolist()
        df['cosmicTissues'] = df['cosmicTissues'].str.split('|').tolist()
        df['cosmicFreqs'] = df['cosmicFreqs'].astype(str).str.split('|').tolist()
        df['cosmicSampleCounts'] = df['cosmicSampleCounts'].astype(str).str.split('|').tolist()

        ## split out the data over multiple rows where a single variant
        ## has multiple cosmic ids or tissue types

        # create an empty list that will hold the new row data
        rows = []

        # create new rows from existing ones depending on number of cosmic ids in the original row
        df.apply(
            lambda row: [
                rows.append(
                    [
                        row['id'],
                        row['cosmicIds'][i],
                        row['cosmicTissues'][i],
                        row['cosmicFreqs'][i],
                        row['cosmicSampleCounts'][i]
                    ]
                )
                for i in range(len(row['cosmicIds']))], axis=1)

        # reassign df by creating a 'new' df from the list of rows
        df = pd.DataFrame(rows, columns=df.columns)

        # convert cosmicFreqs from str to float
        df['cosmicFreqs'] = pd.to_numeric(df['cosmicFreqs'])

        # convert cosmicSampleCounts from str to int
        df['cosmicSampleCounts'] = pd.to_numeric(df['cosmicSampleCounts'])

        # calculate number of samples by multiplying cosmicFreqs and cosmicSampleCounts
        df['cosmicSampleNumber'] = (df['cosmicFreqs']*df['cosmicSampleCounts']).round().astype(int)

        # drop unwanted columns
        df.drop(['cosmicFreqs', 'cosmicSampleCounts'], axis=1, inplace=True)

        ## the aim of the step is to make the variants that share a COSMIC ID but have
        ## multiple, non-haem tissue types all 'look' the same so that we can drop them

        # put a long string into a short variable
        hl = 'Haematopoietic and lymphoid tissue'

        # set cosmicTissues for non-haem variant to NaN
        df.loc[df['cosmicTissues'] != hl, 'cosmicTissues'] = 'Non-haem'

        # set cosmicSampleNumber for non-haem variant to 0
        df.loc[df['cosmicTissues'] != hl, 'cosmicSampleNumber'] = 0

        # drop duplicates
        df.drop_duplicates(inplace=True)

        ## we still have a situation where a variant can be listed twice under a single
        ## cosmic id - once with haem&lymph, and once as non-haem. in this case, we want
        ## to drop the non-haem row, and leave behind only the haem&lymph row

        # count the number of times that each variant is present in the df
        c = Counter(df['id'])

        # create a list of variants that appear more than once
        dups = [n for n in c if c[n] > 1]

        # create a mask for variants in df_new that appear more than once
        mask = df['id'].isin(dups)

        # keep rows that aren't duplicated, or which are duplicated but are in haem&lymph
        df = df.loc[(~mask) | ((mask) & (df['cosmicTissues'] == hl))]

        ## we still have the problem of the variants that have multiple cosmic ids
        ## that are ALL reported in haem&lymph tissue - in this case, we will combine
        ## the data into pipe-separated lists, as in the original alamut batch output

        # sort first by id and second by sample number
        df.sort_values(
            ['id', 'cosmicSampleNumber'],
            ascending = [True, False],
            axis = 0,
            inplace = True
        )

        # re-count the variants in df to look for duplicates
        c = Counter(df['id'])

        # re-assign duplicated variants to dups
        dups = [n for n in c if c[n] > 1]

        ##########################################################################################
        # for some unexplained reason, the next bit of code causes pandas to give a
        # 'SettingWithCopyWarning' error - it is acknowledged by the pandas devs that
        # this warning can be a 'false positive' and doesn't necessarily mean that
        # there's anything wrong with the code.
        #
        # the warning is manually turned off/on again at the start/end of the for loop
        #
        # references:
        # https://github.com/pandas-dev/pandas/issues/10954
        # http://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
        ##########################################################################################

        # loop over the list of duplicates and put the relevant values into pipe-separated lists
        for d in dups:

            pd.options.mode.chained_assignment = None
            df['cosmicIds'].loc[df['id'] == d] = '|'.join(df['cosmicIds'].loc[df['id'] == d].tolist())
            df['cosmicSampleNumber'].loc[df['id'] == d] = '|'.join(df['cosmicSampleNumber'].loc[df['id'] == d].astype(str).tolist())
            pd.options.mode.chained_assignment = 'warn'

        # drop the duplicate rows created in the previous step
        df.drop_duplicates(inplace=True)

        # drop the cosmicTissues column since we're only interested in haem&lymph anyway
        df.drop('cosmicTissues', axis=1, inplace=True)

        # rename some columns
        df.rename(
            columns={
                'cosmicIds': 'cosmic_id',
                'cosmicSampleNumber': 'cosmic_haem'
            },
            inplace=True
        )

        # convert cosmic_haem from int->str to stop it converting to float in the output df
        df['cosmic_haem'] = df['cosmic_haem'].astype(str)

        ## put the 'new' cosmic data back into the original df

        # specify columns to drop from the original df
        cols = ['cosmicIds', 'cosmicTissues', 'cosmicFreqs', 'cosmicSampleCounts']

        # drop the cosmic columns from the original df
        df_whole.drop(cols, axis=1, inplace=True)

        # insert the cosmic columns from the new df
        df_whole = df_whole.merge(df, how='left', left_on='id', right_on='id')

        # change the 'keep' flag to True for variants with a cosmic id
        df_whole.loc[
            (pd.notnull(df_whole['cosmic_id'])) &
            (df_whole['cosmic_haem'] != '0') &
            (df_whole['FREQ'] >= 0.05) &
            (df_whole['DP'] >= 200),
            'keep'
        ] = True

        return df_whole
    
    
    def splice(self, df):
        
        '''Adds splicing annotations based on Alamut Batch data'''
        
        # flag potential splice variants
        df.loc[
            (df['distNearestSS'].between(-3, 3, inclusive=True)) &
            (df['nearestSSChange'] <= -0.2),
            'splice'
        ] = 'splice_effect'

        # create a list of possible consequences
        consequences = [
            'synonymous',
            'missense',
            'stop_gain',
            'in-frame',
            'frameshift',
            'start_loss',
            'stop_loss'
        ]

        # add ',splice_effect' to the consequence column where applicable
        for c in consequences:

            df.loc[
                (df['consequence'] == c) &
                (pd.notnull(df['splice'])),
                'consequence'
            ] = c + ',' + df['splice'].astype(str)

        # for those variants without a 'consequence', take this from the 'splice' column
        df.loc[pd.isnull(df['consequence']), 'consequence'] = df['splice']

        # drop unwanted splice-related columns
        df.drop(['distNearestSS', 'nearestSSChange', 'splice'], axis=1, inplace=True)

        return df
    
    
    def filter_variants(self, df):
        
        '''Filters variants on various criteria'''
        
        print('\nFiltering variants')
        
        # drop variants with no consequence
        df = df.loc[pd.notnull(df['consequence'])]
        
        # drop synonymous variants
        df = df.loc[df['consequence'] != 'synonymous']

        # keep deleterious variants with 5-10% VAF and >500 reads
        consequences = [
            'frameshift',
            'stop_gain',
            'splice_effect',
            'in-frame',
            'start_loss'
        ]

        for c in consequences:

            df.loc[
                (df['consequence'].str.contains(c)) &
                (df['FREQ'] >= 0.05) &
                (df['DP'] >= 200),
                'keep'
            ] = True

        ## for some unexplained reason, the next bit of code causes pandas to give a
        ## 'SettingWithCopyWarning' error - it is acknowledged by the pandas devs that
        ## this warning can be a 'false positive' and doesn't necessarily mean that
        ## there's anything wrong with the code.
        ##
        ## this is solved simply by assigning df as a copy of itself. 
        ##
        ## references:
        ## https://github.com/pandas-dev/pandas/issues/10954
        ## http://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas

        # reassign df as a copy of itself (see above for explanation)
        df = df.copy()

        # keep any other variants that are currently 'false' but which have
        # VAF >= 0.1 and aren't synonymous
        df.loc[
            (df['keep'] == False) &
            (df['FREQ'] >= 0.1) &
            (df['DP'] >= 100),
            'keep'
        ] = True

        # discard variants with 1000g_AF > 0.01
        df.loc[
            df['1000g_AF'] >= 0.01,
            'keep'
        ] = False

        # discard variants with ExAC_AF > 0.01
        df.loc[
            (df['exacAltFreq_all'] > 0.01),
            'keep'
        ] = False

        # keep only the variants where 'keep' is True
        df = df.loc[df['keep'] == True]

        # drop the 'keep' column
        df.drop('keep', axis=1, inplace=True)

        return df
    
    
    def tweak(self, df):
        
        '''Makes various small changes to the output prior to saving'''
        
        ### add the new 'variant' column ###
        
        # convert 'POS' from int to str
        df['POS'] = df['POS'].astype(str)

        # convert e.g. 'chr1' -> '1'
        df['chr'] = df['chr'].str.replace('chr', '')

        # specify the columns with which to create the ID
        cols = ['gene', 'chr', 'POS', 'REF']

        # create an ID for each variant by concatenating gene, chr, pos and ref
        df['variant'] = df[cols].apply(lambda x: '_'.join(x), axis=1)
        
        # add '/ALT' to the end
        df['variant'] = df['variant'] + '/' + df['ALT']

        # specify the columns to drop
        cols = ['chr', 'POS', 'REF', 'ALT', 'id']

        # drop unwanted columns
        df.drop(cols, axis=1, inplace=True)
        
        ### change some columns names ###
        
        # create a dict of column names to change
        cols = {
            'DP': 'depth',
            'FREQ': 'vaf',
            'varType': 'var_type',
            'rsId': 'dbsnp',
            'exacAltFreq_all': 'exac_AF'
        }

        # change column names
        df.rename(columns=cols, inplace=True)

        ### change the order of columns in the output file ###
        
        # specify the order of columns
        cols = [
            'sample_id',
            'gene',
            'variant',
            'vaf',
            'depth',
            'var_type',
            'consequence',
            'cNomen',
            'pNomen',
            'cosmic_id',
            'cosmic_haem',
            'dbsnp',
            '1000g_AF',
            'exac_AF'
        ]

        # update the column order
        df = df[cols]

        ### final tweaks ###
        
        # get around 'SettingWithCopyWarning' for the umpteenth time
        df = df.copy()

        # sort by variant
        df.sort_values('variant', axis=0, inplace=True)

        # reset the index
        df.reset_index(drop=True, inplace=True)

        return df
    

if __name__ == "__main__":
    args = get_args()
    generator = ResultsGenerator(args)
    generator.run()