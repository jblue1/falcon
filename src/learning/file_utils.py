from datetime import date
import os
import shutil
import json
import numpy as np
import pandas as pd
import time
import sys


def make_save_directory(loss_function):
    """
    Creates directory to save losses and weights for each run
    returns - path to the directory
    """
    today = str(date.today())
    run_number = 0
    save_dir = "/Run_" + loss_function + "_" + str(run_number) + "_" + today

    while os.path.exists(save_dir):
        run_number += 1
        save_dir = "/Run_" + loss_function + "_" + str(run_number) + "_" + today

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert os.path.isdir(save_dir)
    return save_dir


def save_network(save_dir, model_path):
    """
    Copy the model.py file into the save directory so weights can easily be loaded
    for future training/evaluations
    save_dir - location to save file
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
    with open(fname, "w") as f:
        json.dump(params_dict, f, indent=4)



def save_losses(save_dir, losses_dict):
    """
    Save loss curves to a txt
    save_dir - location to save file
    losses_dict - dictionary with keys being the name of the
    losses, and the values being a list with the losses
    """
    fname = os.path.join(save_dir, "losses.txt")
    loss_df = pd.DataFrame.from_dict(losses_dict)
    loss_df.to_csv(fname, header="None", index="None", sep=" ")



def main():
    pass


if __name__ == "__main__":
    main()


    

