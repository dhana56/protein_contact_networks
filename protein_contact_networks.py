import pandas as pd
import os
from biopandas.pdb import PandasPdb
from scipy.spatial import distance
import networkx as nx


def pdb_retrive(pdb_id):
    """Function returns the pdb file as a Pandas Dataframe.

    parameters:
    pdb_id: protein PDB ID"""
    if os.path.isfile(pdb_id):
        ppdb = PandasPdb().read_pdb(pdb_id)
    else:
        ppdb = PandasPdb().fetch_pdb(pdb_id)
    return ppdb


def PCN_fun(pdb_file, chain,dir=os.getcwd() , cut_off=7.00, atom='CA', residue_no_diff=0 ):
    """Function used to calculate the protein contact networks. 
    
    parameters:
        pdb_file: pdb_file of the protein. it only accept .pdb format
                  pass the of the file.
        dir: directory_name, to which output get stored.
        chain:    chain  name of PDB ID.
        cut_off:  7.00 angstrom (default)
        atom:     c-alpha atoms of amino acid (default) 
        residue_no_diff: 1, difference between 2 residue atoms 

    Return:
    CSV file, it inlcudes the edges of the protein contact network."""
    
    try:
        #pdb_file handling or retrieve pdb file
        ppdb =  pdb_retrive(pdb_file)
    
        #converting the pdb_file into pandas dataframe.
        residue_df1 = ppdb.df['ATOM'] # residue atoms in the pdb file 
        residue_df =residue_df1.loc[residue_df1['chain_id']==chain.upper() ]
        calpha_df = residue_df[residue_df.atom_name.isin([atom])].reset_index(drop= True)
        
        if len(calpha_df)!=0:

            #Dropping rows the alt_loc column value rather than A- based on the B-factor conformation changes.
            calpha_df_A2 = calpha_df.drop(calpha_df.index[(calpha_df["alt_loc"]=="B")| (calpha_df["alt_loc"]=="C")])
            c_residue_number = calpha_df_A2['residue_number'].values
            
            #extracting the x,y,z  cordiantes of the c-alpha resides.
            values1=list(zip(c_residue_number, calpha_df_A2[['x_coord','y_coord','z_coord']].values))
            
            #keeping the residue conatcts based on the eucleadian cut-off distance of 7 Angstrom.
            ty = [(i[0],j[0]) for i in list(values1) for j in list(values1 )
                    if distance.euclidean(i[1],j[1])<=cut_off and i[0]!=j[0] and abs(i[0]-j[0])!=residue_no_diff]
            ty_1 = set(map(lambda x: tuple(sorted(x)),ty))

            # making the edges from the protein contact work # Is it necessary  
            edge_df = pd.DataFrame(ty_1,columns= ["Node1","Node2"])
            edge_df =edge_df.sort_values(by='Node1').reset_index(drop = True)
            M = nx.from_pandas_edgelist(edge_df[:],'Node1', 'Node2')
            edge_df.to_csv(os.path.join(dir,pdb_file+'.csv'))
            del ty,ty_1,edge_df,ppdb,values1 
            print(f"your files are saved as a {pdb_file}.csv. Go and check it out!")
        
        else:
            print(f"Please check the {pdb_file} has the chain id-> '{chain}' is available or not")

    except Exception as e: 
        print(e)
        print(f"Entry is invalid")
        print(f"please check the PDB id: '{pdb_file}' at https://www.rcsb.org/")

