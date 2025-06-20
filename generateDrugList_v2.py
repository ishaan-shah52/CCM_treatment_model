#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Kaitlyn Wintruba, adapted from Anders Nelson and Ani Chandrabhatla 
# Created: 06/03/2024
# Updated: 06/20/2024 fixed multiple drug target issue and cleaned up code

# Screens DrugBank to find drugs that target genes correpsonding to nodes in the network model

import pandas as pd
from urllib.request import Request, urlopen
import re

# INPUTS: 
all_drugs_CSV = pd.read_csv('alldrugbank.csv') # from DrugBank website
DrugBank_names_CSV= pd.read_csv('drugbankvocabulary.csv') # from DrugBank website
network_model_EXCEL = pd.read_excel('CCMNetIS.xlsx', skiprows=[0]) 
# NOTE: model file needs to have 'gene name' col with multiple genes for a given node separated by ';'

# OUTPUT: DrugsToSimulate.csv

    
#%% IDmatch: match DrugBank drug IDs to the gene names of nodes in the model

def IDmatch(targets, model):

    # Create empty lists to store values
    drugIDs = []; geneSymbols = []; geneNames = []; geneList = []
    
    # Get list of genes in the network model
    model_genes=model['Gene Name'].dropna()
    for gene in model_genes:
        genes = gene.strip().split(';')
        for g in genes:
            geneList.append(g)
            
    # Match gene names to DrugBank drug IDs
    for i in geneList:
        genes_matched = targets[targets['Gene Name'].str.lower() == i.lower()]
        if not genes_matched.empty:
            IDs = genes_matched['Drug IDs'].values[0].split(';')
            for drug_ID in IDs:
                drugIDs.append(drug_ID)
                geneSymbols.append(i)
                geneNames.append(genes_matched['Name'].values[0])
                
    
    # Create dataframe of drug IDs, gene names, and gene symbols
    drug_DF = pd.DataFrame({'DrugID':drugIDs,'Target':geneNames,'Gene':geneSymbols})
    
    # Ensure that the 'DrugID' values in drug_DF have no extra whitespace
    drugIDs_strip = pd.Series([0]*len(drug_DF))
    for i in range(0,len(drug_DF)):
        drugIDs_strip[i] = drug_DF['DrugID'][i].strip()
    drug_DF['DrugID'] = drugIDs_strip
    
    # IDmatch function returns drug_DF
    return drug_DF

drug_DF = IDmatch(all_drugs_CSV,network_model_EXCEL)
print('Drug ID match done!')


#%% webScrapeDrugAction: get drug action info from DrugBank website

def webScrapeDrugAction(drug_DF):
    
    # Create empty list to store values
    actions = []; approved = []
    
    # Collect drug action info for a given drug ID
    for y in range(0, len(drug_DF)):
        
        # Get the drug target and ID from drug_DF
        target = drug_DF.Target[y]
        ID = drug_DF.DrugID[y]; print(ID)
        
        # Scrape all the info on the DrugBank website for a drug ID and compile into string list
        url = ('https://www.drugbank.ca/drugs/'+ID)
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        webUrl = urlopen(req)
        page = webUrl.read().decode('utf-8') # decode converts bytes to string object
        
        # Remove empty rows from scraped page
        striphtml = lambda data: re.split(re.compile(r'<.*?>'), data)
        stripped_page = striphtml(page)
        stripped = list(filter(lambda x: len(x) > 1, stripped_page))
        
        # Find if 'Approved' is indicated on the stripped page
        group = [i for i, x in enumerate(stripped) if x == 'Approved']
        if len(group) > 1: approved.append('Yes')
        else: approved.append('No')
        
        # Find the indices where the target name appears on the stripped page
        indices = [i for i, x in enumerate(stripped) if x == target]
        
        # The word 'Kind' comes after the target name in the targets box on the DrugBank
        # page, so determine where this word follows the target name
        keep = []
        for k in indices:
            if stripped[k+1]=='Kind':
                keep.append(k)
         
        # Check that exactly 1 index was saved to keep
        if(len(keep)!=1):
            print('Warning:', target+' not found on', ID, 'webpage!')
            actions.append('Error')
        
        # Search for drug action following correct target name index
        if(len(keep)==1):
            ant_bool=False
            ag_bool=False
            
            for i in range (keep[0],keep[0]+10):
                word=stripped[i]
                if (word.lower() == 'inhibitor') or (word.lower()=='antagonist') or (word.lower()=='antagonists' or (word.lower() == 'inhibitors')):
                    ant_bool=True
                    break
                else: ant_bool=False   
                    
            for i in range (keep[0],keep[0]+10):   
                word=stripped[i]
                if (word.lower() == 'agonist') or (word.lower() == 'agonists') or (word.lower() == 'inducer') :
                    ag_bool=True
                    break
                else: ag_bool=False 
            
            # Assign drug action values for each drug
            if (ant_bool & ~ag_bool): actions.append('Antagonist')
            elif (~ant_bool & ag_bool): actions.append('Agonist')
            elif (~ant_bool & ~ag_bool): actions.append('Not Available')
            else: actions.append('Conflicting')    
    
    # Add drug action and FDA approval values to drug_DF
    drug_DF['Action']=actions
    drug_DF['FDA Approved']=approved
    
    # webScrapeDrugAction function returns updated drug_DF
    return drug_DF

drug_DF = webScrapeDrugAction(drug_DF)
print('Drug action done!')


#%% IDpull: compile drug info into CSV file

def IDpull(targets, drug_DF, model, DB):
    
    # Create empty lists to store values
    nodeIDs = []; drugnames = []; isagonist = []; agonisttarget = []
    agonisttargetgeneid = []; agonisttargetindex = []; isantagonist = []
    antagonisttarget = []; antagonisttargetgeneid = []; antagonisttargetindex = []

    # Assign node ID values for each gene
    nodeGene = list(model['Gene Name'])
    for y in range(0, len(drug_DF)):
        if drug_DF.Gene[y] in nodeGene:
            i = nodeGene.index(drug_DF.Gene[y])
        else:
            for z in range(0, len(nodeGene)): 
                if type(nodeGene[z]) == float: continue
                elif drug_DF.Gene[y] in nodeGene[z].split(';'): i = z
        nodeIDs.append(list(model['ID'])[i])
      
    # Assign common drug name values for each drug ID
    DB_drugs = list(DB['DrugBank ID']) 
    for y in range(0, len(drug_DF)):
        i = DB_drugs.index(drug_DF.DrugID[y])
        drugnames.append(list(DB['Common name'])[i])
    
    # Add drug name and node name values to drug_DF
    drug_DF['Drug Name'] = drugnames
    drug_DF['Node Name'] = nodeIDs
    
    # Fill in data lists with agonist and antagonist info
    for i in range(0,len(drug_DF)):
        if drug_DF['Action'][i] == 'Agonist':
            isagonist.append('Yes')
            isantagonist.append('No')
            agonisttarget.append(nodeIDs[i])
            agonisttargetgeneid.append(list(drug_DF['Gene'])[i])
            agonisttargetindex.append(str(list(model['ID']).index(nodeIDs[i])+1))
            antagonisttarget.append(''); antagonisttargetgeneid.append(''); antagonisttargetindex.append('')
        elif drug_DF['Action'][i] == 'Antagonist':
            isagonist.append('No')
            isantagonist.append('Yes')
            antagonisttarget.append(nodeIDs[i])
            antagonisttargetgeneid.append(list(drug_DF['Gene'])[i])
            antagonisttargetindex.append(str(list(model['ID']).index(nodeIDs[i])+1))
            agonisttarget.append(''); agonisttargetgeneid.append(''); agonisttargetindex.append('')
        else:
            isagonist.append(''); isantagonist.append(''); agonisttarget.append(''); agonisttargetgeneid.append(''); 
            agonisttargetindex.append(''); antagonisttarget.append(''); antagonisttargetgeneid.append('');
            antagonisttargetindex.append(''); 
        
    # Create Drugs to Simulate Matrix (DSM) as dataframe     
    ds = {'Drug':tuple(drugnames),'IsAgonist':tuple(isagonist),'AgonistTargetGeneID':tuple(agonisttargetgeneid),
          'AgonistTarget':tuple(agonisttarget),'AgonistTargetIndex':tuple(agonisttargetindex),'IsAntagonist': tuple(isantagonist),
          'AntagonistTargetGeneID':tuple(antagonisttargetgeneid),'AntagonistTarget':tuple(antagonisttarget),'AntagonistTargetIndex':tuple(antagonisttargetindex),
          'DrugAction':tuple(['Competitive']*len(drug_DF)), 'DrugApproved':tuple(list(drug_DF['FDA Approved']))}
    DSM = pd.DataFrame(ds)
    # NOTE: Competitive binding is assumed for all drugs and will need to be verified manually
    
    # Remove drugs with no drug action info and sort by drug
    DSM = DSM[DSM['IsAgonist'] != ''].sort_values(by='Drug').reset_index(drop=True)
    
    # Combine multiple rows in DSM for the same drug
    check = lambda series: series.iloc[0] if series.nunique() == 1 else ";".join(series.unique())
    acheck = lambda series: series.iloc[0] if series.nunique() == 1 else 'Yes'
    DSM = DSM.groupby(by='Drug').agg({'Drug':'first',
                                      'IsAgonist':acheck,
                                      'AgonistTargetGeneID':check,
                                      'AgonistTarget':check,
                                      'AgonistTargetIndex':check,
                                      'IsAntagonist':acheck,
                                      'AntagonistTargetGeneID':check,
                                      'AntagonistTarget':check,
                                      'AntagonistTargetIndex':check,
                                      'DrugApproved':'first',
                                     'DrugAction':'first'})
    
    # Clean up DSM by removing extra ';'
    clean = lambda x: x.replace(';;', ';').strip('; ').lstrip(';')
    DSM = DSM.applymap(lambda x: clean(x) if isinstance(x, str) else x)
 
    # Write DSM to CSV file
    DSM.to_csv('DrugsToSimulate.csv', index = False, header = True)
    
    # IDpull function returns drugs to simulate matrix dataframe
    return DSM
    
DSM = IDpull(all_drugs_CSV, drug_DF, network_model_EXCEL, DrugBank_names_CSV)
print('DrugsToSimulate.csv created!')