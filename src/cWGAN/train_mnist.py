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
    label = np.zeros(10)
    label[number] = 1
    return label


def main():
    today = str(date.today())
    run_number = 0
    save_dir = '../../models/cWGAN/Run_MNIST_' + today + '_' + str(run_number) 
    
    while (os.path.exists(save_dir)):
        run_number += 1
        save_dir = '../../models/cWGAN/Run_MNIST' + today + '_' + str(run_number) 

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert(os.path.isdir(save_dir))


    (train_images, train_labels), (_, _) = keras.datasets.mnist.load_data()
    train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
    train_images = (train_images - 127.5) / 127.5 # Normalize the images to [-1, 1]

    one_hot_labels = np.zeros((train_labels.shape[0], 10))
    for i in range(len(train_labels)):
        one_hot_labels[i, :] = one_hot_encode(train_labels[i])

    dataset = tf.data.Dataset.from_tensor_slices((one_hot_labels, train_images))

    print(dataset)
    cwgan = model.cWGAN_mnist(train.NUM_CRITIC_ITERS, train.BATCH_SIZE,
            noise_dims=100)
    fname = os.path.join(save_dir, 'modelsummary.txt')
    with open(fname, 'w') as f:
        with redirect_stdout(f):
            cwgan.print_network()


    start = time.time()
    losses, critic_losses = train.train(cwgan, dataset, train.NUM_EPOCHS, train_labels.shape[0], save_dir)
    print("Time for {} epochs: {}".format(train.NUM_EPOCHS, time.time() - start))
    filename = save_dir + '/losses.txt'
    critic_filename = save_dir + '/critic_losses.txt'
    loss_df = pd.DataFrame({"Wass Est": pd.Series(losses)})
    loss_df.to_csv(filename, header='None', index='None', sep=' ')
    critic_loss_df = pd.DataFrame({"Critic_loss": pd.Series(critic_losses)})
    critic_loss_df.to_csv(critic_filename, header='None', index='None', sep=' ')


if __name__ == '__main__':
    main()
