import os
from datetime import datetime
from ete3 import NCBITaxa
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)
ncbi = NCBITaxa()
sqldb = os.path.join(os.environ['HOME'], '.etetoolkit', 'taxa.sqlite')  # path to ete sql db
db_modification_time = datetime.fromtimestamp(os.path.getctime(sqldb))
database_age_days = abs((db_modification_time-datetime.now()).days)
if database_age_days > 30:
    ncbi.update_taxonomy_database()
    comment = 'taxa.sqlite is more than a month old, updating database. More details in logs/benchmark'
else:
    comment = 'taxa.sqlite is up to date'
file = open(snakemake.output[0], 'w')
file.write(comment)
file.close()