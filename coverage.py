import argparse
import os
import pandas as pd


def get_args():
    
    desc = '''
    1) Generates a coverage statement for the sample, 
    2) Generates a per amplicon coverage file for the sample, 
    3) Optionally generates a hotspot coverage file for the sample.
    '''
    
    parser = argparse.ArgumentParser(
        description = desc
    )
    
    parser.add_argument(
        '--pileup',
        help = 'Path to pileup file',
        metavar = '<input.pileup>',
        dest = 'pileup',
        required = True
    )
    
    parser.add_argument(
        '--bed',
        help = 'Path to BED file',
        metavar = '<input.bed>',
        dest = 'bedfile',
        required = True,
    )
    
    parser.add_argument(
        '--hotspots',
        help = 'Path to hotspots file',
        metavar = '<hotspots.txt>',
        dest = 'hotspots',
        required = False,
        default = False
    )
    
    parser.add_argument(
        '--outdir',
        help = 'Path to output directory',
        metavar = '<outdir>',
        dest = 'outdir',
        required = True
    )
    
    parser.add_argument(
        '--srsf2',
        help = 'Use this flag to include SRSF2 in the coverage statement',
        dest = 'srsf2',
        action = 'store_true',
        default = False,
        required = False
    )
    
    args = parser.parse_args()
    
    return args


class CoverageGenerator(object):
    
    def __init__(self, args):
        
        # grab the command line arguments
        self.pileup = args.pileup
        self.bedfile = args.bedfile
        self.hotspots = args.hotspots
        self.outdir = args.outdir
        self.srsf2 = args.srsf2
        
        # grab the sample name from the pileup file name
        self.sample = os.path.basename(self.pileup).replace('.pileup', '')
        
        return None
    
    
    def run(self):
        
        self.check_files()
        self.setup()
        
        df_pile = self.load_pileup()
        
        self.generate_coverage_statement(df_pile)

        self.generate_amplicon_coverage(df_pile)
        
        if self.hotspots:
            self.generate_hotspot_coverage(df_pile)
        
        return None
    
    
    def check_files(self):
        
        '''
        Checks that the specified input files exist
        Checks that the specified input files have correct suffixes
        Checks that the specified output directory exists
        '''
        
        # check that input files exist
        files = {
            'pileup file': self.pileup,
            'BED file': self.bedfile
        }
        
        if self.hotspots:
            files['Hotspots file'] = self.hotspots
        
        for k, v in files.iteritems():
            if os.path.isfile(v):
                print('\nOK: %s exists: %s' % (k, v))
            else:
                sys.exit('\nERROR: %s does not exist: %s' % (k, v))
                
        # check that the specified input files have correct suffixes
        suffixes = {
            self.pileup: '.pileup',
            self.bedfile: '.bed'
        }
        
        if self.hotspots:
            suffixes[self.hotspots] = '.txt'
        
        for k, v in suffixes.iteritems():
            if not k.endswith(v):
                sys.exit('\nERROR: %s does not have a %s suffix' % (k, v))

        # check that the specified output directory exists
        if os.path.isdir(self.outdir):
            print('\nOK: Output directory exists: %s' % self.outdir)
        else:
            sys.exit('\nERROR: Output directory does not exist: %s' % self.outdir)
        
        return None
    
    
    def setup(self):
        
        '''
        Make subdirectories for coverage statements, amplicon coverage and hotspot coverage
        '''
        
        subdirs = ['coverage_statements', 'amplicon_coverage', 'hotspot_coverage']
        
        os.chdir(self.outdir)
        
        for subdir in subdirs:
            
            if not os.path.isdir(subdir):
            
                os.mkdir(subdir)
            
        return None
        
        
    def load_pileup(self):

        '''
        Loads a subset of columns from the pileup file into a Pandas dataframe
        before returning the df
        '''
        
        print('\nLoading pileup file')
        
        df_pile = pd.read_table(
            self.pileup,
            names = ['chr', 'pos', 'base', 'count', 'foo', 'bar'],
            usecols = ['chr', 'pos', 'count']
        )

        return df_pile
    
    
    def load_bedfile(self):

        '''
        Loads a subset of columns from the bed file into a Pandas dataframe and
        drops any duplicate amplicons before returning the df
        '''
        
        print('\nLoading BED file')
        
        df_bed = pd.read_table(
            self.bedfile,
            names = ['chr', 'start', 'end', 'strand', 'amplicon'],
            usecols = ['chr', 'start', 'end', 'amplicon']
        )

        # correct the zero-based numbering of the start position
        df_bed['start'] = df_bed['start'] + 1

        # use drop duplicates to drop any amplicons with identical start/end
        df_bed.drop_duplicates(subset=['start', 'end'], inplace=True)

        # reset the index
        df_bed = df_bed.reset_index(drop=True)

        return df_bed
    

    def generate_coverage_statement(self, df_pile):

        '''
        Generates a coverage statement for the sample
        '''
        
        meancount = int(df_pile['count'].mean())
        numbases = float(len(df_pile))
        over100 = float(len(df_pile.loc[df_pile['count'] > 100]))
        over100pct = ((over100 / numbases) * 100)

        pt1 = 'For sample %s mean coverage is %d and %.1f%% of bases have at least 100X coverage. ' % (self.sample, meancount, over100pct)
        pt2 = 'The hotspot region of SRSF2 is poorly covered and will be resequenced. '
        pt3 = 'See datafiles for full mutational profile. '
        
        # if the --srsf2 flag has been used
        if self.srsf2:
        
            # attempt to create a dataframe from the srsf2 hotspot
            df_srsf2 = df_pile.loc[
                (df_pile['chr'] == 'chr17') &
                (df_pile['pos'] == 74732959)
            ]

            # if this yields an empty df, the hotspot must not be covered
            if df_srsf2.empty:
                # therefore set srsf2 to zero
                srsf2 = 0
            else:
                # otherwise set srsf2 to the value of 'count'
                srsf2 = df_srsf2['count'].iloc[0]

            # having set a value to srsf2, if it's <20
            if srsf2 < 20:
                # include the srsf2 section in the coverage statement
                coverage_statement = pt1 + pt2 + pt3
            else:
                # otherwise, leave it out
                coverage_statement = pt1 + pt3
                
        else:
            coverage_statement = pt1 + pt3
        
        outfile = os.path.join(
            self.outdir,
            'coverage_statements',
            self.sample + '_coverageStatement.txt'
        )
        
        with open(outfile, 'w') as f:
            f.write(coverage_statement)
        
        print('\nCoverage statement saved to %s' % outfile)

        return None
    
    
    def generate_amplicon_coverage(self, df_pile):
        
        '''
        Generates amplicon coverage file for the sample
        '''
        
        # load the bed file into a dataframe
        df_bed = self.load_bedfile()
        
        # get a list of chromosomes featured in the bed file
        chromlist = self.get_chromlist(df_bed)
        
        # create an empty dataframe to be filled with each iteration over chromlist
        df_all = pd.DataFrame(columns = ['amplicon', 'avgcount'])
        
        print('\nGenerating amplicon coverage\n')
        
        for chrom in chromlist:
    
            print('Processing chromosome %s' % chrom.replace('chr',''))

            # subset df_bed based on chrom
            df_tmp = df_bed.loc[df_bed['chr'] == chrom]

            intervals = self.get_intervals(df_tmp)

            pairlist = self.get_pairlist(intervals)

            df = self.get_region_counts(pairlist, df_pile, df_tmp)

            # for those regions only covered by a single amplicon, the coverage for that amplicon
            # is equal to the coverage for that region. create df_out from these regions in df
            df_out = df.loc[df['numamps'] == 1]
            df_out = df_out[['amplicon', 'avgcount']]

            # initially the amplicon names are comma-separated strings;
            # convert these into lists so they can be iterated over
            df['amplicon'] = df['amplicon'].str.split(',')

            df_zero = self.get_df_zero(df)

            # append df_zero to df_out, drop dups, reset index
            df_out = df_out.append(df_zero, ignore_index=True)
            df_out.drop_duplicates(inplace=True)
            df_out = df_out.reset_index(drop=True)

            # keep only the rows in df where numamps != 1
            df = df.loc[df['numamps'] != 1]

            # keep only the rows in df where avgcount != 0
            df = df.loc[df['avgcount'] != 0]

            while not df.empty:

                df = self.check_completed_amplicons(df, df_out)

                df_out = self.calculate_amplicon_coverage(df, df_out)

            df_all = df_all.append(df_out, ignore_index=True)

        # for some reason avgcount becomes a float, so change back to int
        df_all['avgcount'] = df_all['avgcount'].astype(int)

        # merge df_bed with df_all
        df_all = pd.merge(df_bed, df_all, how='inner', on='amplicon')
        
        # ensure that any negative values are changed to zero
        df_all.loc[df_all['avgcount'] < 0, 'avgcount'] = 0
        
        # set the output file name
        outfile = os.path.join(
            self.outdir,
            'amplicon_coverage',
            self.sample + '_ampliconCoverage.txt'
        )

        # save as tab-separated file
        df_all.to_csv(outfile, sep='\t', header=True, index=False)

        print('\nAmplicon coverage file saved to %s' % outfile)
        
        return None
    
    
    def get_chromlist(self, df_bed):
        
        '''
        Gets a sorted list of chromosomes represented in the BED file
        '''
        
        chromlist = [i.replace('chr', '') for i in list(set(df_bed['chr'].tolist()))]
        chromlist = sorted([24 if i=='X' else int(i) for i in chromlist])
        chromlist = ['chrX' if i==24 else 'chr' + str(i) for i in chromlist]
        
        return chromlist
    
    
    def get_intervals(self, df):

        '''Returns a list of amplicon start and end points'''

        return sorted(list(set(df['start'].tolist() + [i+1 for i in df['end'].tolist()])))
    
    
    def get_pairlist(self, intervals):

        pairlist = []

        for count, item in enumerate(intervals):

            if count < (len(intervals)-1):

                pairlist.append([intervals[count], intervals[count+1]-1])

        return pairlist
    
    
    def get_region_counts(self, pairlist, df_pile, df_tmp):

        for pair in pairlist:

            # calculate average coverage over that region and append to 'pair'
            avgcount = df_pile['count'].loc[df_pile['pos'].isin(range(pair[0], pair[1]+1))].mean()
            pair.append(avgcount)

            # get the start and end values for the pair
            startA = pair[0]
            endA = pair[1]

            # start an empty list to hold the names of overlapping amplicons
            amplist = []

            # loop over df_bed checking if that amplicon overlaps with the region of interest
            for i in range(len(df_tmp)):

                # grab the amplicon name
                amp = df_tmp['amplicon'].iloc[i]

                # set the start and end values for the amplicon
                startB = df_tmp['start'].iloc[i]
                endB = df_tmp['end'].iloc[i]

                # if there is overlap, add the amplicon name to the list
                if (startA <= endB) and (endA >= startB):
                    amplist.append(amp)

            pair.append(','.join(amplist))
            pair.append(len(amplist))

        # put pairlist into a df
        df = pd.DataFrame(
            pairlist,
            columns = ['start', 'end', 'avgcount', 'amplicon', 'numamps']
        )

        # fill NaN values as 0 and convert float -> int
        df['avgcount'].fillna(value=0, inplace=True)
        df['avgcount'] = df['avgcount'].astype(int)

        # drop any regions that aren't covered by any amplicons
        df = df.loc[df['numamps'] > 0]

        return df
    
    
    def get_df_zero(self, df):

        # create a flat list of amplicon names from what is initially a list of lists
        # when you select df['amplicon'] where df['avgcount'] == 0
        flat_list = list(
            set(
                [item for sublist in df['amplicon'].loc[df['avgcount'] == 0].tolist() for item in sublist]
            )
        )

        # create df_zero from flat_list
        df_zero = pd.DataFrame(flat_list, columns = ['amplicon'])

        # add a new column 'avgcount' for which the value for each row is zero
        df_zero['avgcount'] = 0

        return df_zero
    
    
    def check_completed_amplicons(self, df, df_out):

        '''
        for each row in df, this checks whether ALL of the amplicons
        corresponding to that row are currently represented in df_out.
        it creates a list of True/False values and puts these into a new
        column in df named 'drop', before dropping those rows from df
        '''

        droplist = []

        for i in range(len(df)):

            if all(item in df_out['amplicon'].tolist() for item in df['amplicon'].iloc[i]):

                droplist.append(True)

            else:

                droplist.append(False)

        df['drop'] = droplist

        df = df.loc[df['drop'] == False]

        df = df.reset_index(drop=True)

        return df
    
    
    def calculate_amplicon_coverage(self, df, df_out):

        for i in range(len(df)):

            # create a list of amplicons for that row that are represented in df_out
            matching = [amp for amp in df['amplicon'].iloc[i] if amp in df_out['amplicon'].tolist()]

            # create a list of amplicons for that row that are NOT represented in df_out
            notmatching = [amp for amp in df['amplicon'].iloc[i] if amp not in df_out['amplicon'].tolist()]

            # if the number of matching amplicons is equal to numamps-1 then
            # that means there is only a single amplicon from that row for which
            # the avgcount value has not yet been calculated.
            if len(matching) == df['numamps'].iloc[i] - 1:

                # df_tmp is created from df_out based on the matching amplicon values
                df_tmp = df_out.loc[df_out['amplicon'].isin(matching)]

                # the result for the single 'notmatching' amplicon for that row
                # is equal to the 'avgcount' value for that region minus the sum
                # of 'avgcount' values for the other 'matching amplicons'
                result = df['avgcount'].iloc[i] - df_tmp['avgcount'].sum()

                # df_result is created to hold the 'nonmatching' amplicon name
                # and avgcount value and this is then appended to df_tmp
                df_result = pd.DataFrame([[notmatching[0], result]], columns = list(df_tmp))

                # append df_result to df_tmp
                df_tmp = df_tmp.append(df_result)

                # append df_tmp to df_out
                df_out = df_out.append(df_tmp)

                # drop duplicates and reset index
                df_out.drop_duplicates(inplace=True)
                df_out = df_out.reset_index(drop=True)

        return df_out
    
    
    def generate_hotspot_coverage(self, df_pile):
        
        '''
        Generates hotspot coverage file for the sample
        '''
        
        # load hotspots
        df_hotspots = pd.read_table(self.hotspots, usecols=['chr', 'POS', 'gene'])
        
        # rename 'POS' to 'pos'
        df_hotspots.rename(columns={'POS': 'pos'}, inplace=True)

        # drop dups and reset index
        df_hotspots.drop_duplicates(inplace=True)
        df_hotspots = df_hotspots.reset_index(drop=True)
        
        # carry out a left join of df_hotspots and df_pile
        df = pd.merge(df_hotspots, df_pile, on=['chr', 'pos'], how='left')

        # fill NaN values with 0
        df['count'].fillna(value=0, inplace=True)

        # ensure df['count'] is an int
        df['count'] = df['count'].astype(int)
        
        outfile = os.path.join(
            self.outdir,
            'hotspot_coverage',
            self.sample + '_hotspotCoverage.txt'
        )

        df.to_csv(outfile, sep='\t', header=True, index=False)

        print('\nHotspot coverage file saved to %s' % outfile)


if __name__ == "__main__":
    args = get_args()
    covgen = CoverageGenerator(args)
    covgen.run()