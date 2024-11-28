import Bio.PDB.PDBParser as PDBParser
import numpy as np
import matplotlib.pyplot as plt
import requests
import os

def download_pdb_file(pdb_id, save_path="."):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_file_path = os.path.join(save_path, f"{pdb_id}.pdb")
        with open(pdb_file_path, "wb") as file:
            file.write(response.content)
        print(f"PDB file {pdb_id} downloaded successfully.")
        return pdb_file_path
    else:
        raise Exception(f"Failed to download PDB file {pdb_id} (status code: {response.status_code}).")

def calculate_distance(atom1, atom2):
    diff_vector = atom1.coord - atom2.coord
    return np.sqrt(np.sum(diff_vector ** 2))

def contact_map(structure, distance_threshold=8.0):
    ca_atoms = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    
    num_atoms = len(ca_atoms)
    contact_matrix = np.zeros((num_atoms, num_atoms))
    
    for i, atom1 in enumerate(ca_atoms):
        for j, atom2 in enumerate(ca_atoms):
            if i != j:
                distance = calculate_distance(atom1, atom2)
                if distance <= distance_threshold:
                    contact_matrix[i, j] = 1
                    
    return contact_matrix

def plot_contact_map(contact_matrix):
    plt.figure(figsize=(10, 10))
    plt.imshow(contact_matrix, cmap='Greys', origin='upper')
    plt.colorbar(label='Contact')
    plt.title('Contact Map')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.show()

pdb_id = "4YWO"
pdb_file_path = download_pdb_file(pdb_id)

parser = PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdb_file_path)

contact_matrix = contact_map(structure, distance_threshold=8.0)

plot_contact_map(contact_matrix)