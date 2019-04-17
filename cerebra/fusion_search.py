""" searches STAR-fusion output files (.tsv) for a given fusion of 
	interest. implements standard two-gene run mode, as well as single
	gene query.  """

import pandas as pd
import os
import sys
import click


def search_func_any(row, GOI):
	""" search for ANY partner to a given gene """
	cellFile = row['name']
	cellName = cellFile.split('_')[0] + '_' + cellFile.split('_')[1]

	PATH = cwd + 'fusions/' + cellFile
	
	curr_fusions = pd.read_csv(PATH, sep='\t')

	try:
		fusionsList = list(curr_fusions['#FusionName'])
	except KeyError:
		print('error caught, in %s' % cellName)
		outputRow = pd.DataFrame([[cellName, 'ERROR']])
		return outputRow 

	outputRow = pd.DataFrame([[cellName, 0]])
		
	for item in fusionsList:
		item0 = item.split('--')[0]
		item1 = item.split('--')[1]
		if GOI == item0 or GOI == item1:
			print(cellName)
			print(item)
			outputRow = pd.DataFrame([[cellName, 1]])
			return outputRow

	return outputRow



def search_func(row, FOI):
	""" standard search for a given two-gene fusion. """
	cellFile = row['name']
	cellName = cellFile.split('_')[0] + '_' + cellFile.split('_')[1]

	PATH = cwd + 'fusions/' + cellFile
	
	curr_fusions = pd.read_csv(PATH, sep='\t')
	
	try:
		fusionsList = list(curr_fusions['#FusionName'])
	except KeyError:
		print('error caught, in %s' % cellName)
		outputRow = pd.DataFrame([[cellName, 'ERROR']])
		return outputRow 

	if str(FOI) in fusionsList:
		outputRow = pd.DataFrame([[cellName, 1]])
		print(cellName)
	else:
		outputRow = pd.DataFrame([[cellName, 0]])

	return outputRow



""" get cmdline input """
@click.command()
@click.option('--test', default = False)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)
 


def fusion_search(test, wrkdir):
	""" searches STAR-fusion output files for a particular fusion of interest """
	global colNames
	global cwd 

	test_bool = test
	cwd = wrkdir

	print(' ')
	print(' ')

	cmd = 'sudo mkdir -p ' + cwd + 'fusionOut'
	cmd1 = 'sudo chmod -R 777 ' + cwd + 'fusionOut'
	os.system(cmd)
	os.system(cmd1)

	colNames = ['cellName', 'fusionPresent_bool']

	cellFiles = os.listdir(cwd + 'fusions')
	cellFiles_df = pd.DataFrame(data=cellFiles, columns=['name']) # need to convert to df before apply call

	print('running...')

	with open(cwd + '../fusionsList.csv', 'r') as f:
		fusionsList = f.readlines()
		fusionsList = [x.strip() for x in fusionsList]

		for currFus in fusionsList:
			outFileStr = cwd + 'fusionOut/' + currFus + '.query.out.csv'

			queryStrSplit = currFus.split('--')
			queryStrRev = queryStrSplit[1] + '--' + queryStrSplit[0]

			outputRows = cellFiles_df.apply(search_func, axis=1, args=(currFus,))
			outputRows_list = list(outputRows) # for some reason need to convert to a list before concatting
			outputDF = pd.concat(outputRows_list, ignore_index=True)

			outputRows_rev = cellFiles_df.apply(search_func, axis=1, args=(queryStrRev,))
			outputRows_rev_list = list(outputRows_rev) # for some reason need to convert to a list before concatting
			outputDF_rev = pd.concat(outputRows_rev_list, ignore_index=True)

			outputDF.columns = colNames
			outputDF_rev.columns = colNames

			outputDF_comb = pd.DataFrame(columns=colNames)
			outputDF_comb['cellName'] = outputDF['cellName']
			outputDF_comb['fusionPresent_bool'] = outputDF['fusionPresent_bool'] + outputDF_rev['fusionPresent_bool']
			outputDF_comb.to_csv(outFileStr, index=False)

	print('done!')
