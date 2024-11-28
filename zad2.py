import os
import requests
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt

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

def load_structure(pdb_id, save_path="."):
    pdb_file_path = download_pdb_file(pdb_id, save_path)
    structure = PDB.PDBParser(QUIET=True).get_structure(pdb_id, pdb_file_path)
    return structure

def calculate_phi_psi(structure):
    model = structure[0]  
    phi_psi_angles = []
    
    for chain in model:
        polypeptides = PDB.PPBuilder().build_peptides(chain)
        for poly in polypeptides:
            phi_psi = poly.get_phi_psi_list()
            phi_psi_angles.extend(phi_psi)
    
    phi_psi_angles = [(phi, psi) for phi, psi in phi_psi if phi is not None and psi is not None]
    return phi_psi_angles

def plot_ramachandran(phi_psi_angles):
    phi_angles, psi_angles = zip(*phi_psi_angles)  # Rozpakowujemy listę
    phi_angles = np.degrees(phi_angles)
    psi_angles = np.degrees(psi_angles)

    plt.figure(figsize=(8, 8))
    plt.scatter(phi_angles, psi_angles, c='blue', alpha=0.5, s=10)
    plt.title("Wykres Ramachandrana")
    plt.xlabel("Phi (°)")
    plt.ylabel("Psi (°)")
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.grid(color='lightgray', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.show()

# Główna część skryptu
if __name__ == "__main__":
    pdb_id = "4YWO"
    save_path = "."  # Katalog do zapisu pliku PDB
    structure = load_structure(pdb_id, save_path)
    phi_psi_angles = calculate_phi_psi(structure)
    plot_ramachandran(phi_psi_angles)