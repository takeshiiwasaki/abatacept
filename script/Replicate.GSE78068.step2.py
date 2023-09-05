import sys
import click
import numpy as np
import pandas as pd
import os
os.environ['R_HOME'] = 'C:\Program Files\R\R-4.1.0'
import rpy2.robjects as ro
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from numpy.testing import _private
from statsmodels.stats.multitest import multipletests

#Preparation
p=pd.read_table("C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/pheno.xls")
e=pd.read_table("C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/expr.xls")

p["ID"]=p.index
p.index=range(len(p))
for i in range(len(p)):
    if "Non" in p.loc[i,"title"]:
        p.loc[i,"Target"]="zero"
    else:
        p.loc[i,"Target"]="one"
p=p[["ID","Target"]]
p.to_csv("C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/design.xls",sep="\t",index=False)

##Gene Annotatation for expression file
a = pd.read_table("C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/annot.GSE78068.xls")
a = a[["agilent_wholegenome_4x44k_v2","external_gene_name"]]
a = a.rename(columns={'agilent_wholegenome_4x44k_v2':'ID'})
e["ID"] = e.index
m=pd.merge(e,a,on="ID",how="outer").dropna().drop_duplicates().reset_index(drop=True)
m=m.rename(columns={'external_gene_name':'Symbol'})
m=m.drop('ID', axis=1)
e_final=pd.DataFrame()
for g in set(m["Symbol"].tolist()):
    tmp=pd.DataFrame({g:m[(m["Symbol"]==g)].mean()}).T
    e_final=pd.concat([e_final,tmp])

e_final.to_csv("C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/expr_target_annotate.xls",sep="\t",index=True)

#Read input file
data = pd.read_table('C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/expr_target_annotate.xls')
data=data.rename(columns={'Unnamed: 0': 'ID'})
data=data.set_index("ID")
design = pd.read_table('C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/design.xls') #replace with your own design file

#Import R libraries
base = importr('base')
stats = importr('stats')
limma = importr('limma')
writexl = importr('writexl')
    
# Convert data and design pandas dataframes to R dataframes
with localconverter(ro.default_converter + pandas2ri.converter):
    #r_data = ro.conversion.py2ri(data)
    r_data=pandas2ri.py2rpy(data)
    #r_design = ro.conversion.py2ri(design)
    r_design= pandas2ri.py2rpy(design)
    # Use the genes index column from data as a R String Vector
    genes = ro.StrVector(
        [
            str(index)
            #added tovalues to convert to numpy array
            for index in data.index.tolist()
            #for index in data.index.tolist()
        ]
    )

# Create a model matrix using design's Target column using the R formula "~0 + f" to get all the unique factors as columns
f = base.factor(r_design.rx2('Target'), levels=base.unique(r_design.rx2('Target')))
form = Formula('~0 + f')
form.environment['f'] = f
r_design = stats.model_matrix(form)
r_design.colnames = base.levels(f)

# Fit the data to the design using lmFit from limma
fit = limma.lmFit(r_data, r_design)
# Make a contrasts matrix with the 1st and the last unique values
contrast_matrix = limma.makeContrasts(f"{r_design.colnames[0]}-{r_design.colnames[-1]}", levels=r_design)

# Fit the contrasts matrix to the lmFit data & calculate the bayesian fit
fit2 = limma.contrasts_fit(fit, contrast_matrix)
fit2 = limma.eBayes(fit2)

# topTreat the bayesian fit using the contrasts and add the genelist
r_output = limma.topTreat(fit2, coef=1, genelist=genes, number=np.Inf)
writexl.write_xlsx(r_output, "C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/limma_output.xlsx")
pd.read_excel("C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/limma_output.xlsx").to_csv("C:/Users/tiwasaki/Documents/RA_RNA.ABT/GSE78068/limma_output.xls",sep="\t",index=False)
