import sys
import cGAN
import file_utils
import data_utils
import FCNN
import cWGAN
import tensorflow as tf
import os
import time
import numpy as np


def train_fcnn(params):
    print("Training fcnn")
    params_dict = file_utils.get_FCNN_hyperparams(params)
    trainer = FCNN.Trainer(params_dict)
    history = trainer.train()
    trainer.save_losses(history)
    trainer.save_model()
    trainer.save_params(params_dict)


def train_cGAN(params):
    print("Training cGAN")
    params_dict = file_utils.get_cGAN_hyperparams(params)
    trainer = cGAN.Trainer(params_dict)
    trainer.train()
    trainer.save_losses()
    trainer.save_model()
    trainer.save_params(params_dict)


def train_cWGAN(params):
    print("Training cWGAN")
    params_dict = file_utils.get_cWGAN_hyperparams(params)
    trainer = None
    if params_dict["data_path"] == "MNIST":
        trainer = cWGAN.MNISTTrainer(params_dict)
    else:
        trainer = cWGAN.Trainer(params_dict)
    trainer.train()
    trainer.save_losses()
    trainer.save_model()
    trainer.save_params(params_dict)


def train(model, params):
    if model == "FCNN":
        train_fcnn(params)

    elif model == "cGAN":
        train_cGAN(params)

    elif model == "cWGAN":
        train_cWGAN(params)

    else:
        usage("Invalid model chosen")


def usage(message):
    print(message)
    print("usage: python train.py model param_file")
    print("    model - Model to train. One of FCNN, cGAN, cWGAN")
    print("    param_file - path to json file with model parameters")
    sys.exit(1)


def main():
    if len(sys.argv) != 3:
        usage("Incorrect number of arguments given")
    model = sys.argv[1]
    params = sys.argv[2]
    train(model, params)


if __name__ == "__main__":
    main()