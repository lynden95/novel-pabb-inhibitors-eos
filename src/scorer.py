#import all libraries 
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdmolfiles import MolToPDBFile
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from standardiser import standardise



#define pathnames
DATAPATH = "../data/"
RESULTSPATH = "../results/"
SOURCEPATH = "../src"

#define standard arguments/parameters

#smiles prep - standardise and filter
smiles_csv = os.path.join(DATAPATH, "smiles", "smiles.csv")
filtered_std_smiles = os.path.join(DATAPATH, "smiles", "filtered_std_smiles.csv") 

#ligand prep
sdf_folder = os.path.join(DATAPATH, "smiles")
pH = 7.4

#docking
receptor = os.path.join(DATAPATH, "protein", "pabb_model1.pdbqt") 
ligands = os.path.join(DATAPATH, "smiles", "merged.sdf")  
logfile = os.path.join(RESULTSPATH, "outputs", "log.txt") 
docking_output = os.path.join(RESULTSPATH, "outputs", "outputs.sdf") 
smina = os.path.join(SOURCEPATH, "smina.static") 

#analysis
affinities_id_smiles = os.path.join(RESULTSPATH, "outputs", "affinity_id_smiles.csv")



#define functions

def prep_smiles(smiles_csv):

    df=pd.read_csv(smiles_csv) 

    mols = [Chem.MolFromSmiles(smi) for smi in df["SMILES"].tolist()]

    std_mols = []

    for mol in mols:
        if mol is not None:
            try:
                std_mol = standardise.run(mol)
            except:
                std_mol = np.nan
        else:
            std_mol = np.nan
        std_mols += [std_mol]

    std_smiles = []

    for std_mol in std_mols:
        if std_mol is not None:
            try: 
                std_smi = Chem.MolToSmiles(std_mol)
            except:
                std_smi=np.nan
        else:
            std_smi = np.nan
        std_smiles += [std_smi]

    df["ST_SMILES"] = std_smiles

    df=df[df["ST_SMILES"].notna()]

    df.drop(columns = ["SMILES"], inplace=True)

    mols = [Chem.MolFromSmiles(smi) for smi in df["ST_SMILES"].tolist()]

    mol_wts = []

    for mol in mols:
        mol_wt = Chem.Descriptors.MolWt(mol)
        mol_wts += [mol_wt]

    df["mol_wt"] = mol_wts

    df_filtered = df[df['mol_wt'] >= 300]
    df_filtered = df_filtered[df_filtered['mol_wt'] <= 800]
    df_filtered.drop(columns = ["mol_wt"], inplace=True)
    df_filtered.to_csv(filtered_std_smiles, index=False)



def prepare_ligands_sdf(filtered_std_smiles, pH, header_len=1, output_dir = sdf_folder, delim=',') -> list:

    out_sdfs = [] 
    
    print(filtered_std_smiles)
    
    with open(filtered_std_smiles, 'r') as csv: 
            
        for entry in csv.readlines()[header_len:]:
            
            ID, ST_SMILES = entry.split(delim)[:2]            
           
            # Convert smiles str to 3D coordinates
            mol = Chem.MolFromSmiles(ST_SMILES)
            mol = Chem.AddHs(mol)
            ps = AllChem.ETKDGv2()
            ps.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, ps)
                        
            if type(result) is int:
                if result == -1:
                    continue
           
            # Ouput coords to pdb
            pdb_name = f"{output_dir}/{ID}.pdb"
            MolToPDBFile(mol, pdb_name)
       
            # Protonate according to pH, convert to .sdf
            sdf_name = f"{output_dir}/{ID}.sdf"
            ! obabel {pdb_name} -pH {pH} -O {sdf_name} 
          
            os.remove(pdb_name) # removes the .pdb files after obabel protonates and converts to .sdf
                       
            out_sdfs.append(sdf_name)

            os.remove(filtered_std_smiles)

    return out_sdfs



def merge_sdfs(out_sdfs, merged_sdf):

    mols = []
    for s in out_sdfs:
        suppl = Chem.SDMolSupplier(s)
        for  mol in suppl:
            if mol is None:
                continue
        else:
            mols += [mol]
        os.remove(s)

    with Chem.SDWriter(merged_sdf) as w:
        for mol in mols:
            if mol is None:
                continue
            else:
                w.write(mol)



def prepare_and_merge_ligands(filtered_std_smiles,  pH, header_len=1, output_dir='', delim=','):
    ligands = prepare_ligands_sdf(filtered_std_smiles,  pH, header_len, output_dir, delim)
    merge_sdfs(ligands, output_dir+"/merged.sdf") 



def docking_cmd():
    smina  + " -r " + receptor + " -l " + ligands + " -o " + docking_output + " --log " + logfile + " --seed 42 " + " --center_x 74 " + " --center_y 44 " + " --center_z 57 " + " --size_x 22 " + " --size_y 22 " + " --size_z 24 " + " --exhaustiveness 8 " + " --num_modes 1 " + " --addH off " + " --scoring vinardo "


def affinities_to_smiles(docking_output):

    df = PandasTools.LoadSDF(docking_output, embedProps=True, molColName=None, smilesName='SMILES')
    df['ID'] = df['ID'].map(lambda x: x.lstrip('../data/smiles/').rstrip('.pdb'))
    df.to_csv(affinities_id_smiles, index=False)


def run():
    prep_smiles(smiles_csv)
    prepare_and_merge_ligands(filtered_std_smiles, pH, output_dir=sdf_folder)
    os.remove(filtered_std_smiles)
    os.system(docking_cmd)
    os.remove(ligands)
    affinities_to_smiles(docking_output)


    
    