""" creates fusions x cell dataframe, as well as fusions x patient
    needed for input to makeSummaryTable script """

import numpy as np
import pandas as pd


def readFunc(fus):
	""" reads the fusions csv and looks for non-zero entries """
	fName = cwd + 'fusionOut/' + fus + '.query.out.csv'
	f = pd.read_csv(fName)
	toKeep = f['fusionPresent_bool'] == 1
	f = f[toKeep]
	f = f.reset_index(drop=True)
	
	return f



""" get cmdline input """
@click.command()
@click.option('--test', default = False)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)



def fusions_x_cell(test, wrkdir):
	""" create a cell-wise fusion counts table. 
		needed for creating validationTable"""

	global cwd

	cwd = wrkdir 

	with open(cwd + 'fusionsList.csv', 'r') as f:
		fusionsList = f.readlines()
		fusionsList = [x.strip() for x in fusionsList]

		fusionsDF = pd.DataFrame('Nan', index=np.arange(50), columns=fusionsList)

		for currFus in fusionsList:
    		currName = currFus.split('--')[0] + '_' + currFus.split('--')[1]
    		df = readFunc(currFus)
    		fusionsDF[currFus] = pd.Series(df['cellName'])

		fusionsDF.to_csv(cwd + 'fusion_dataframe.csv', index=False)
