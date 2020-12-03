#!/usr/bin/env python3

# Import libraries
import sys
import numpy as np
import pandas as pd

# Get file name
in_file = sys.argv[1]

# Read file as data frame
df = pd.read_csv(in_file, header=None)

# Take single columns
singles = df.loc[:, 177:180]

# Remove columns with empty entries
df_new = df.dropna(axis = 1)

# Reset column indices
df_new.columns = np.arange(len(df_new.columns))

# Replace entries
pd.set_option('mode.chained_assignment', None)
df_new.loc[1, 1:22] = "ref"
df_new.loc[1, 23:45] = "test"

# Bind for final DF
out = pd.concat([df_new, singles], axis = 1, ignore_index = True)

# Overwrite file
out.to_csv(in_file, header = False, index = False)
