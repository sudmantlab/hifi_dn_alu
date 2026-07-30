### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen'].unique()
specimens_by_group = sample_table.groupby('group')['specimen'].unique().apply(list).to_dict()

# workaround right now for not being able to predefine chromosomes/contigs for parallelization of DeepVariant
chrs = ['chr' + str(n) for n in np.arange(1, 23).tolist()+['X', 'Y']]
