import pkgutil
import pandas as pd
from io import BytesIO

def get_phenotype_table():
    b = BytesIO(pkgutil.get_data(__name__, 'phenotype-table.csv'))
    return pd.read_csv(b)
