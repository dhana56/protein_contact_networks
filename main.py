from protein_contact_networks import PCN_fun

if __name__=="__main__":

    pdb_id = input("Enter the PDB id here -> ")
    chain_id = input("Enter the chain id here -> ")
    PCN_fun(pdb_file=pdb_id,chain= chain_id)
    