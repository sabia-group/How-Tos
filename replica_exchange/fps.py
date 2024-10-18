# This script is to do the farthest point sampling and select the representative structures for the training set.
# Use as python fps.py -i inputfile -n nos
# its ase read input and ouput, so can be used for any format of ase readable file

# For this script, you need to have the following packages installed:
# librascal, scikit-learn, ase, numpy, matplotlib, chemiscope, skmatter


import numpy as np
import argparse

from ase.io import read, write

#import chemiscope
from rascal.representations import SphericalInvariants as SOAP
from skmatter.preprocessing import StandardFlexibleScaler
from sklearn.decomposition import PCA
from skmatter.feature_selection import FPS
#from matplotlib import pyplot as plt

# set the plotting style
#plt.rcParams['axes.linewidth'] = 4
#plt.rcParams['xtick.major.size'] = 10
#plt.rcParams['ytick.major.size'] = 10
#plt.rcParams['xtick.major.width'] = 2
#plt.rcParams['ytick.major.width'] = 2
#plt.rcParams['xtick.direction'] = 'in'
#plt.rcParams['ytick.direction'] = 'in'
#plt.rcParams['font.family'] = 'sans-serif'  # Set font family to sans-serif
#plt.rcParams['font.sans-serif'] = 'Open Sans'  # Set font to Open Sans
#plt.rcParams['font.size'] = 22
#plt.rcParams['mathtext.fontset'] = 'custom'
#plt.rcParams['mathtext.it'] = 'STIXGeneral:italic'
#plt.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
#plt.rcParams["figure.figsize"] = (8, 8)
#plt.rcParams["legend.framealpha"] = 1
#plt.rcParams["legend.fancybox"] = False


def main(inputfile, nos):
    structure_file = inputfile
    frames = read(structure_file, index=':', format='extxyz')

    print('Number of frames: ', len(frames), flush=True)
    print('Number of atoms/frame: ', len(frames[0]), flush=True)

    SOAP_HYPERS = {
        "interaction_cutoff": 3.5,
        "max_radial": 6,
        "max_angular": 6,
        "gaussian_sigma_constant": 0.4,
        "cutoff_smooth_width": 0.5,
        "gaussian_sigma_type": "Constant",
    }
    print("SOAP hyperparameters\n: ", SOAP_HYPERS, flush=True)

    numbers = list(sorted(set([int(n) for frame in frames for n in frame.numbers])))

    # initialize SOAP

    soap = SOAP(
        global_species=numbers, 
        expansion_by_species_method='user defined',
        **SOAP_HYPERS            
    )

    X = None 
    print("computing SOAP features...", flush=True)
    for i, frame in enumerate(frames):
        # normalize cell for librascal input
        if np.linalg.norm(frame.cell) < 1e-16:
            extend = 1.5 * (np.max(frame.positions.flatten()) - np.min(frame.positions.flatten()))
            frame.cell = [extend, extend, extend]
            frame.pbc = True
        frame.wrap(eps=1e-16)

        x = soap.transform(frame).get_features(soap).mean(axis=0) # here it takes mean over atoms in the frame
        if X is None:
            X = np.zeros((len(frames), x.shape[-1]))
        X[i] = x

    print(f"SOAP features shape: {X.shape}", flush=True)

    # np.save('full-featurization.npy', X)
    X_full = X.copy()
    n_FPS = nos # number of structures to select 
    print(f"Performing farthest point sampling -------\nPlease wait a while", flush=True)
    struct_idx = FPS(n_to_select=n_FPS, initialize = 123).fit(X.T).selected_idx_
    X_fps = X[struct_idx]

    print(f"FPS selected indices: {struct_idx.shape}")
    print(f"Original: {X.shape} ---> FPS: {X_fps.shape}")

    frames_fps = [frames[i] for i in struct_idx]
    print(f"fps selected structures are in {inputfile.split('.')[0]}-fps.extxyz file", flush=True)
    write(f"{inputfile.split('.')[0]}-fps.extxyz", frames_fps, format='extxyz')

    frames_not_fps = [frames[i] for i in range(len(frames)) if i not in struct_idx]
    print(f"rest number of structures are in {inputfile.split('.')[0]}-rest.extxyz file", flush=True)
    print("use rest file for pool", flush=True)
    write(f"{inputfile.split('.')[0]}-rest.extxyz", frames_not_fps, format='extxyz')

    #energy_full = [frame.info['energy'] for frame in frames]
    ## PCA on FPS selections
    # X_full = StandardFlexibleScaler(column_wise=False).fit_transform(X_full)
    # T = PCA(n_components=2).fit_transform(X_full)
    #plt.figure(figsize=(10, 8))
    #scatter = plt.scatter(T[:,0], T[:,1], c=np.array(energy_full).reshape(-1, 1), s=100, cmap='viridis', alpha=0.7, label = 'Full Dataset')
    #cbar = plt.colorbar(scatter)
    #cbar.ax.tick_params(labelsize=20)
    #cbar.set_label('Energy (eV)', fontsize=20)
    ##plt.scatter(X[:,0], X[:,1], s=100, alpha=0.7, c='grey', label = 'Full Dataset')
    #plt.scatter(T[struct_idx,0], T[struct_idx,1], c='red', s=100 ,label = 'FPS Selected', facecolors='none', edgecolors='black')
    #plt.title('PCA on SOAP features')
    #plt.xlabel('PCA 1')
    #plt.ylabel('PCA 2')
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(f"{inputfile.split('.')[0]}-fps.pdf", dpi=300, bbox_inches='tight')
    #plt.close()

    ## energy plot
    #plt.figure(figsize=(8, 8))
    #plt.plot(energy_full,'go' ,label='Full Dataset')
    #plt.plot(struct_idx, np.array(energy_full)[struct_idx], 'ro', label='FPS Selected')
    #plt.xlabel('Structure index')
    #plt.ylabel('Energy (eV)')
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(f"{inputfile.split('.')[0]}-energy.pdf", dpi=300, bbox_inches='tight')
    #plt.close()

    ## forces plot
    #forces_full = [frame.arrays['forces'] for frame in frames]
    #forces_fps = [frame.arrays['forces'] for frame in frames_fps]
    #plt.figure(figsize=(8, 8))
    #plt.plot([np.linalg.norm(f) for f in forces_full],'go' ,label='Full Dataset')
    #plt.plot(struct_idx, [np.linalg.norm(f) for f in forces_fps], 'ro', label='FPS Selected')
    #plt.xlabel('Structure index')
    #plt.ylabel('Force (eV/A)')
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(f"{inputfile.split('.')[0]}-forces.pdf", dpi=300, bbox_inches='tight')
    #plt.close()


if __name__ =="__main__":
    parser = argparse.ArgumentParser( description = """This script is to do the farthest point sampling on the SOAP features and select the representative structures for the training set.""")

    parser.add_argument('-i', type = str, help='input data file')
    parser.add_argument('-n', type = int, help='number of fps structures to select')


    args = parser.parse_args()
    inputfile = args.i
    nos = args.n

    main(inputfile, nos)
