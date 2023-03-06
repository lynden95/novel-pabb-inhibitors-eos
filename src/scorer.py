#import all libraries 
import os
import pandas as pd
import numpy as np
import rdkit
import standardiser
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdmolfiles import MolToPDBFile
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from standardiser import standardise # do I need to have the standardiser notebook sent from GT in my notebooks?

#define pathnames
DATAPATH = "../data/"
RESULTSPATH = "../results/"
SOURCEPATH = "../src"

#define standard arguments/parameters

#standardising
csv_file = os.path.join(DATAPATH, "smiles", "smiles.csv")
std_csv_file = os.path.join(DATAPATH, "smiles", "std_smiles.csv")

#filtering
csv_200_descriptors = os.path.join(DATAPATH, "smiles", "csv_200_desc.csv")
smiles_plus_desc = os.path.join(DATAPATH, "smiles", "smiles_descriptors.csv") 
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
affinities = os.path.join(RESULTSPATH, "outputs", 'affinities.csv')
affinities_header = os.path.join(RESULTSPATH, "outputs", "affinities_header.csv")
final_smiles = os.path.join(DATAPATH, "smiles", "final_merged_smiles.csv")
affinities_and_smiles = os.path.join(RESULTSPATH, "outputs", "affinities_with_stdsmiles.csv")



#define functions
def smiles_standardiser(csv_file, std_csv_file):

    df=pd.read_csv(csv_file) 

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

    df.to_csv(std_csv_file, index=False)

    #os.remove(csv_file) can remove the original file if we don't need to keep it!



def mol_descriptor(std_csv_file):
    
    mols = [Chem.MolFromSmiles(i) for i in std_csv_file] 
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()
    
    Mol_descriptors =[]
    for mol in mols:
        # add hydrogens to molecules
        mol=Chem.AddHs(mol)
        # Calculate all 200 descriptors for each molecule
        descriptors = calc.CalcDescriptors(mol)
        Mol_descriptors.append(descriptors)
    return Mol_descriptors,desc_names 



def filterer(std_csv_file, filtered_std_smiles):
    df = pd.read_csv(std_csv_file)

    Mol_descriptors,desc_names = mol_descriptor(df1['ST_SMILES'])

    df_200_descriptors = pd.DataFrame(Mol_descriptors,columns=desc_names)

    df_200_descriptors.to_csv(csv_200_descriptors, index=False)

    cols = [5]
    df_molwt = df_200_descriptors[df_200_descriptors.columns[cols]]

    merged = pd.concat([df, df_molwt], axis="columns")

    merged.to_csv(smiles_plus_desc, index=False)
    os.remove(csv_200_descriptors)

    df = pd.read_csv(smiles_plus_desc)

    # Filter all rows for which the smiles mw is under 100
    df_filtered = df[df['MolWt'] >= 300]
    df_filtered = df_filtered[df_filtered['MolWt'] <= 800]

    df_filtered.drop(columns = ["MolWt"], inplace=True)

    df_filtered.to_csv(filtered_std_smiles, index=False)

    os.remove(smiles_plus_desc)
    #os.remove(std_csv_file)



def prepare_ligands_sdf(filtered_std_smiles, pH, header_len=1, output_dir = sdf_folder, delim=',') -> list:



def merge_sdfs(out_sdfs, merged_sdf):



def prepare_and_merge_ligands(filtered_std_smiles,  pH, header_len=1, output_dir='', delim=','):
    ligands = prepare_ligands_sdf(filtered_std_smiles,  pH, header_len, output_dir, delim)
    merge_sdfs(ligands, output_dir+"/merged.sdf")



def docking_cmd():
    smina  + " -r " + receptor + " -l " + ligands + " -o " + docking_output + " --log " + logfile + " --seed 42 " + " --center_x 74 " + " --center_y 44 " + " --center_z 57 " + " --size_x 22 " + " --size_y 22 " + " --size_z 24 " + " --exhaustiveness 8 " + " --num_modes 1 " + " --addH off " + " --scoring vinardo "



def affinities_to_smiles(logfile, final_smiles, affinities_and_smiles):

    with open(logfile, "r") as f:
        i = 0;
        for line in f:
            if '-+' in line:
                nextline = next(f)
                i = i + 1

                nextlinearray  = nextline.split()                       #splitting the first row in different values
                bind_aff = nextlinearray[1]                             #getting the binding affinity of first pose

                with open(affinities, "a") as myfile:
                    print(bind_aff, end='\n', file=myfile)
        

    df = pd.read_csv(affinities, names = ['AFFINITY'])

    headerList = ['AFFINITY']
    df.to_csv(affinities_header, header=headerList, index=False)
    df = pd.read_csv(affinities_header)

    os.remove(affinities)

    df1 = pd.read_csv(final_smiles)
    df2 = pd.read_csv(affinities_header)

    merged = pd.concat([df1, df2], axis="columns")
    # merged.drop(columns = ["ID"], inplace=True) -- this removes the ID line, so output is smiles and affininty only... could keep.

    merged.to_csv(affinities_and_smiles, index=False)

    os.remove(affinities_header)



def run():
    os.system(smiles_standardiser(csv_file, std_csv_file))
    os.system(mol_descriptor(std_csv_file))
    os.system(filterer(std_csv_file, filtered_std_smiles))
    os.system(prepare_and_merge_ligands(filtered_std_smiles, pH, output_dir=sdf_folder))
    os.system(docking_cmd())
    os.system(affinities_to_smiles(logfile, final_smiles, affinities_and_smiles))


    
    