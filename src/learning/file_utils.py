from datetime import date
import os
import shutil
import json
import pandas as pd


def make_save_directory(model):
    """Create directory to save losses and weights for each run

    Args:
        model (string): Kind of model being trained (FCNN, cGAN, etc.)

    Returns:
        path-like: Path to directory where all of the training information for a run
        is saved.
    """

    today = str(date.today())
    run_number = 0
    models_dir = "../../models/" + model
    save_dir = models_dir + "/Run_" + today + "_" + str(run_number)

    while os.path.exists(save_dir):
        run_number += 1
        save_dir = models_dir + "/Run_" + today + "_" + str(run_number)

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert os.path.isdir(save_dir)
    return save_dir


def save_network(save_dir, model_path):
    """Copy the file containing the model class into the save directory so weights can
    easily be loaded for future training/evaluation.

    Args:
        save_dir (path-like): Path to the save directory
        model_path (path-like): Path to .py file with the model class
    """

    shutil.copy(model_path, save_dir)


def save_params(save_dir, params_dict):
    """Save the model and training hyperparameters to a json file

    Args:
        save_dir (path-like): Path to save director
        params_dict (dict): Dictionary with model/training parameters
    """

    fname = os.path.join(save_dir, "params.json")
    with open(fname, "w") as f:
        json.dump(params_dict, f, indent=4)


def save_losses(save_dir, losses_dict, prefix=""):
    """Save losses to a txt file

    Args:
        save_dir (path-like): Path to save directory
        losses_dict (dict): Dictionary containing the model losses
        prefix (str, optional): String to append to the file name. Defaults to "".
    """

    name = prefix + "losses.txt"
    fname = os.path.join(save_dir, name)
    loss_df = pd.DataFrame.from_dict(losses_dict)
    loss_df.to_csv(fname, header="None", index="None", sep=" ")


def get_cWGAN_hyperparams(params_file_path):
    """Get hyperparameters from a json file

    Args:
        params_file_path (path-like): Path to json file

    Returns:
        dict: Dictionary with hyperparameters.
    """
    with open(params_file_path, "r") as f:
        params_dict = json.load(f)
    if "cWGAN" in params_dict.keys():
        return params_dict["cWGAN"]
    else:
        return params_dict


def get_FCNN_hyperparams(params_file_path):
    """Get hyperparameters from a json file

    Args:
        params_file_path (path-like): Path to json file

    Returns:
        dict: Dictionary with hyperparameters.
    """
    with open(params_file_path, "r") as f:
        params_dict = json.load(f)
    if "FCNN" in params_dict.keys():
        return params_dict["FCNN"]
    else:
        return params_dict


def get_classifier_hyperparams(params_file_path):
    """Get hyperparameters from a json file

    Args:
        params_file_path (path-like): Path to json file

    Returns:
        dict: Dictionary with hyperparameters.
    """
    with open(params_file_path, "r") as f:
        params_dict = json.load(f)
    if "classifier" in params_dict.keys():
        return params_dict["classifier"]
    else:
        return params_dict


def get_cGAN_hyperparams(params_file_path):
    """Get hyperparameters from a json file

    Args:
        params_file_path (path-like): Path to json file

    Returns:
        dict: Dictionary with hyperparameters.
    """
    with open(params_file_path, "r") as f:
        params_dict = json.load(f)
    if "cGAN" in params_dict.keys():
        return params_dict["cGAN"]
    else:
        return params_dict


def main():
    pass


if __name__ == "__main__":
    main()
