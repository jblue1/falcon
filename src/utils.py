from datetime import date
import os
import shutil
import json
import numpy as np


def make_save_directory():
    """
    Creates directory to save losses and weights for each run
    returns - path to the directory
    """
    today = str(date.today())
    run_number = 0
    save_dir = '../../models/cGAN/Run_' + today + '_' + str(run_number) 
    
    while (os.path.exists(save_dir)):
        run_number += 1
        save_dir = '../../models/cGAN/Run_' + today + '_' + str(run_number) 

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert(os.path.isdir(save_dir))
    return save_dir


def save_network(save_dir, model_path):
    """
    Copy the model.py file into the save directory so weights can easily be loaded
    for future training/evaluations
    save_dir - location to save files
    model_path - path to file containing the model
    """
    shutil.copy(model_path, save_dir)


def save_params(save_dir, params_dict):
    """
    Save the network and training hyper paramters to a json file.
    save_dir - location to save file
    params_dict - dictonary with the hyperparameters
    """
    fname = os.path.join(save_dir, "params.json")
    with open(fname, 'w') as f:
        json.dump(params_dict, f, indent=4)


def load_jet_data(data_path, batch_size):
    """
    Load and normalize the parton data.
    data_path - path to txt file with jet 4-momenta
    returns - tuple of numpy arrays (parton_data, reco_data)
    """
    data = np.loadtxt(data_path, skiprows=2)
    partonPtMax = np.max(data[:, 0], axis=0)
    partonPtMin = np.min(data[:, 0], axis=0)
    partonMean = np.mean(data[:, 1:3], axis=0)
    partonStd = np.std(data[:, 1:3], axis=0)
    partonEMax = np.max(data[:, 3], axis=0)
    partonEMin = np.min(data[:, 3], axis=0)
    
    pfPtMax = np.max(data[:, 4], axis=0)
    pfPtMin = np.min(data[:, 4], axis=0)
    pfMean = np.mean(data[:, 5:7], axis=0)
    pfStd = np.std(data[:, 5:7], axis=0)
    pfEMax = np.max(data[:, 7], axis=0)
    pfEMin = np.min(data[:, 7], axis=0)

    data[:, 0] = (data[:, 0] - partonPtMin)/partonPtMax
    data[:, 1:3] = (data[:, 1:3] - partonMean)/partonStd
    data[:, 3] = (data[:, 3] - partonEMin)/partonEMax
    data[:, 4] = (data[:, 4] - pfPtMin)/pfPtMax
    data[:, 5:7] = (data[:, 5:7] - pfMean)/pfStd
    data[:, 7] = (data[:, 7] - pfEMin)/pfEMax

    np.random.shuffle(data)
    parton_data = data[:, :4]
    reco_data = data[:, 4:]
    return (parton_data, reco_data)

