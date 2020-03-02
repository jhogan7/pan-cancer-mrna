import cptac
import pandas as pd 


def shift(df):
  minVal = (df.min()).min()
  df = df - (minVal -1)
  return df


def get_matched(data):
  patient_list = data.index.values.tolist()
  normal = []
  cancerous = []
  for ID in patient_list:
    if ID.endswith('N'):
      normal.append(ID[:-2])
    else:
      cancerous.append(ID)
  matched = [value for value in normal if value in cancerous]
  return matched

def get_ratio(gene, pairs, mrna, protein):
  cancer_ratio = []
  normal_ratio = []


  for patient in pairs:
    normal = patient + ".N"
    prot_cancer = protein.loc[patient][gene]
    prot_normal = protein.loc[normal][gene]
    mrna_cancer = mrna.loc[patient][gene]
    mrna_normal = mrna.loc[normal][gene]
    cancer_ratio.append(mrna_cancer/prot_cancer)
    normal_ratio.append(mrna_normal/prot_normal)

  return cancer_ratio, normal_ratio

def process_cancer(mrna, protein):
  mrna = shift(mrna)
  protein = shift(protein)
  matched_patients = [value for value in get_matched(protein) if value in get_matched(mrna)]
  gene_list = [value for value in list(protein.columns) if value in list(mrna.columns)] 

  control_data = {}
  cancer_data = {}
  for gene in gene_list: 
    cancer, control = get_ratio(gene, matched_patients, mrna, protein)
    control_data[gene] = control
    cancer_data[gene] = cancer

  cancer_ratios = pd.DataFrame(cancer_data, index=matched_patients)
  control_ratios = pd.DataFrame(control_data, index=matched_patients)

  return cancer_ratios, control_ratios


cptac.download(dataset='Endometrial')
en = cptac.Endometrial()
en_protein = en.get_proteomics()
en_mrna = en.get_transcriptomics()
en_cancer, en_control = process_cancer(en_mrna, en_protein)
en_cancer.to_csv('en_cancer.csv')
en_control.to_csv('en_control.csv')
