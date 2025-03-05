import pandas as pd
import sys
path = sys.argv[1]
df = pd.read_csv(path)
print(df.loc[df['values']<= float(sys.argv[2]), 'relative.size'].max())

