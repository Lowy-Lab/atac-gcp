import pandas as pd
import re
from google.cloud import storage
import sys

client = storage.Client()
sample = sys.argv[2]
exclude = sys.argv[3:]
bucket_name = sys.argv[1].lstrip('gs://').split('/')[0]
bucket = client.get_bucket(bucket_name)
prefix = sys.argv[1].lstrip('gs://').split('/')[1]
prefix = prefix + '/outs/per_sample_outs'
blobs = bucket.list_blobs(prefix=prefix)
x=0
complexities = []
for blob in blobs:
    if re.search(r'library_complexity.csv$', blob.name):
        sample_name = blob.name.split('/')[-2]
        if sample_name in exclude:
            continue
        df=pd.read_csv('gs://' + bucket_name + '/' + blob.name)
        complexities.append(float(df.loc[df['relative.size'] == 1, "values"].iloc[0]))
sample_df = pd.read_csv(f"gs://{bucket_name}/{prefix}/{sample}/library_complexity.csv")
complexity = sample_df.loc[sample_df['values']<=min(complexities), 'relative.size'].max()
print(complexity)
