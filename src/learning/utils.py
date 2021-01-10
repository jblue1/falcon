from datetime import date
import os
import shutil
import json
import numpy as np
import pandas as pd
import time
import sys
import tensorflow as tf


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

    data[:, 0] = (data[:, 0] - partonPtMin) / partonPtMax
    data[:, 1:3] = (data[:, 1:3] - partonMean) / partonStd
    data[:, 3] = (data[:, 3] - partonEMin) / partonEMax
    data[:, 4] = (data[:, 4] - pfPtMin) / pfPtMax
    data[:, 5:7] = (data[:, 5:7] - pfMean) / pfStd
    data[:, 7] = (data[:, 7] - pfEMin) / pfEMax

    np.random.shuffle(data)
    parton_data = data[:, :4]
    reco_data = data[:, 4:]
    return (parton_data, reco_data)


def one_hot_encode(number):
    """
    Return one-hot-encoded vector to serve as a label for MNIST data, for
    example 
    one_hot_encode(4)
    would return
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    number - digit (must be 0 through 9)
    returns - one hot encoded vector
    """
    if (number < 0 or number > 9):
        sys.exit(1)
    label = np.zeros(10)
    label[number] = 1
    return label


def load_mnist_data():
    """
    Load and normalize MNIST data
    returns - tuple of numpy arrays (labels, images)
    """
    (train_images, train_labels), (_, _) = tf.keras.datasets.mnist.load_data()
    train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
    train_images = (train_images - 127.5) / 127.5 # Normalize the images to [-1, 1]

    one_hot_labels = np.zeros((train_labels.shape[0], 10))
    for i in range(len(train_labels)):
        one_hot_labels[i, :] = one_hot_encode(train_labels[i])

    return (train_labels, one_hot_labels, train_images)


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


def train_cGAN(
    dataset,
    model,
    lr,
    batch_size,
    num_critic,
    epochs,    
    use_wasserstein=False,
    save_run=True,
):
    """
    Train an conditional GAN using either a standard or Wasserstein loss
    """
    params_dict = {"lr": lr, "batch_size": batch_size, "num_critic": num_critic}
    if save_run:
        save_dir = ""
        if use_wasserstein:
            params_dict["loss"] = "Wasserstein"
            save_dir = make_save_directory("Wasserstein")
        else:
            params_dict["loss"] = "Normal"
            save_dir = make_save_directory("Normal")

        save_params(save_dir, params_dict)
        checkpoint_dir = save_dir + '/training_checkpoints'


    num_training_examples = len(dataset[0])
    print("Number of training examples: {}".format(num_training_examples))
    for epoch in range(epochs):
        start = time.time()

        indices = np.random.choice(num_training_examples, batch_size, replace=False)
        print("Indices: {}".format(indices))
        labels = dataset[0][indices]
        one_hot_vecs = dataset[1][indices]
        print("Labels")
        print(labels)
        print(one_hot_vecs)


def main():
    dataset = load_mnist_data()
    train_cGAN(dataset, 'hi', 0.1, 5, 1, 1, save_run=False)


if __name__ == "__main__":
    main()


    

