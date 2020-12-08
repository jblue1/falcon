import os
from datetime import date
import time
from contextlib import redirect_stdout
import tensorflow as tf
from tensorflow import keras
import numpy as np
import model
import train
import pandas as pd


def one_hot_encode(number):
    """
    Return 10-d vector of one hot encoded number
    """
    label = np.zeros(10)
    label[number] = 1
    return label


def load_data():
    """
    Load and normalize MNIST data
    returns - tf.data.Dataset object
    returns - number of training examples
    """
    (train_images, train_labels), (_, _) = keras.datasets.mnist.load_data()
    train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
    train_images = (train_images - 127.5) / 127.5 # Normalize the images to [-1, 1]

    one_hot_labels = np.zeros((train_labels.shape[0], 10))
    for i in range(len(train_labels)):
        one_hot_labels[i, :] = one_hot_encode(train_labels[i])

    dataset = tf.data.Dataset.from_tensor_slices((one_hot_labels, train_images))
    return dataset, train_labels.shape[0]


def print_network(save_dir, net):
    """
    Print network layers to text file in save directory
    save_dir - location to save file
    net - cWGAN object 
    """

    fname = os.path.join(save_dir, 'modelsummary.txt')
    with open(fname, 'w') as f:
        with redirect_stdout(f):
            net.print_network()


def main():
    
    save_dir = train.make_save_directory()

    dataset, num_data_points = load_data()
    
    cwgan = model.cWGAN_mnist(train.NUM_CRITIC_ITERS, train.BATCH_SIZE,
            noise_dims=100)
    
    print_network(save_dir, cwgan)
    start = time.time()
    generator_losses, critic_losses = train.train(cwgan, dataset,
            train.NUM_EPOCHS, num_data_points, save_dir)
    print("Time for {} epochs: {}".format(train.NUM_EPOCHS, time.time() - start))
    train.save_losses(save_dir, generator_losses, critic_losses)
    

if __name__ == '__main__':
    main()
