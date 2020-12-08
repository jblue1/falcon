import os
import tensorflow as tf
from tensorflow import keras
import pandas as pd
import numpy as np
import model
from datetime import date
import time
from contextlib import redirect_stdout

NUM_CRITIC_ITERS = 50
BATCH_SIZE = 64
NUM_EPOCHS = 2000


def train(model, dataset, epochs, num_data_points, save_dir):
    '''
    Train model and save weights every 5 epochs. Each training loop requires NUM_CRITIC_ITERS batches of data to train the
    critic, one batch to train the generator, and then one batch to evaluate
    one model performance.
    model - cWGAN object
    dataset - tf.data.Dataset object
    epochs - numer of epochs to train for
    num_data_points - number of training examples in the dataset
    save_dir - location for saving weights
    '''
    dataset = dataset.shuffle(num_data_points)
    dataset = dataset.batch(BATCH_SIZE)
    losses = []
    critic_losses = []
    checkpoint_dir = save_dir + '/training_checkpoints'
    for epoch in range(epochs):
        batches = []
        count = 0
        for batch in dataset:
            count += 1
            batches.append(batch)
            if count % (NUM_CRITIC_ITERS + 2) == 0:
                wass_est, new_critic_losses = model.train_step(batches)
                print("Loss: {}".format(wass_est))
                losses.append(wass_est)
                for loss in new_critic_losses:
                    critic_losses.append(loss)
                batches.clear()
            if count % ((NUM_CRITIC_ITERS + 2)*100) == 0:
                print("Completed {} batches".format(count/(NUM_CRITIC_ITERS+2)))

        if (epoch + 1) % 20 == 0:
            gen_filename = os.path.join(checkpoint_dir, 'gen_' + str(epoch + 1))
            critic_filename = os.path.join(checkpoint_dir, 'critic_' +
                    str(epoch + 1))
            print(gen_filename)
            model.generator.save_weights(gen_filename)
            model.critic.save_weights(critic_filename)

    return losses, critic_losses


def make_save_directory():
    """
    Creates directory to save losses and weights for each run
    returns - path to the directory
    """
    today = str(date.today())
    run_number = 0
    save_dir = '../../models/cWGAN/Run_' + today + '_' + str(run_number) 
    
    while (os.path.exists(save_dir)):
        run_number += 1
        save_dir = '../../models/cWGAN/Run_' + today + '_' + str(run_number) 

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert(os.path.isdir(save_dir))
    return save_dir


def load_data():
    """
    Loads and normalizes data 
    returns - tf.data.Dataset object
    returns - number of training examples
    """
    data = np.loadtxt('../../data/processed/matchedJets.txt', skiprows=2)
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
    num_data_points = len(data)
    np.random.shuffle(data)
    trainParton = data[:, :4]
    trainPf = data[:, 4:]
    dataset = tf.data.Dataset.from_tensor_slices((trainParton,
            trainPf))
    return dataset, num_data_points


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


def save_losses(save_dir, generator_losses, critic_losses):
    """
    Save losses to text file
    save_dir - location to save file
    generator_losses - list of generator losses
    critic_losses - list of critic losses. Note that this
    list should be NUM_CRITIC_ITERS times the length of the
    list of generator losses
    """
    critic_filename = save_dir + '/critic_losses.txt'
    generator_filename = save_dir + '/generator_losses.txt'
    generator_loss_df = pd.DataFrame({"Wass Est": pd.Series(generator_losses)})
    critic_loss_df = pd.DataFrame({"Wass Est": pd.Series(critic_losses)})
    generator_loss_df.to_csv(generator_filename, header='None', index='None', sep=' ')
    critic_loss_df.to_csv(critic_filename, header='None', index='None', sep=' ')


def main():
    save_dir = make_save_directory()
    dataset, num_data_points = load_data()
    
    cwgan = model.cWGAN(NUM_CRITIC_ITERS, BATCH_SIZE)
    print_network(save_dir, cwgan)

    start = time.time()
    generator_losses, critic_losses = train(cwgan, dataset, NUM_EPOCHS, num_data_points, save_dir)
    print("Time for {} epochs: {}".format(NUM_EPOCHS, time.time() - start))
    save_losses(save_dir, generator_losses, critic_losses)
    

if __name__ == "__main__":
    main()

