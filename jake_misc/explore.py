import cptac
from plotnine import *

if __name__ == '__main__':
    import cptac
    datasets=['colon','endometrial','ovarian','cCRcC']
    password_protected_datasets = [ "brca", "luad", "gbm", "hnscc"]
    #[cptac.download(x) for x in datasets]
    cptac.download('endometrial')
    en = cptac.Endometrial()
    print(en.get_cancer_type())
    pt_hist = en.get_medical_history()
    prot = en.get_proteomics()
    phprot = en.get_phosphoproteomics_gene()
    #gene1_filter = acetyl.columns.get_level_values("Name").str.startswith("AA")
    #  Select all columns where the gene starts with "AA". This will grab every column where the key "Name" starts with AA
    #colon.reduce_multiindex(df=phospho, flatten=True)
    #phosphoproteomics has 2 layers
    print()