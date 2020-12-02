import tensorflow as tf
from tensorflow import keras
import pandas as pd
import numpy as np
import model
from datetime import date
import time
from contextlib import redirect_stdout
import os

NUM_CRITIC_ITERS = 50
BATCH_SIZE = 64
NUM_EPOCHS = 10

def train(model, dataset, epochs, num_data_points, save_dir):
    '''
    Each training loop requires NUM_CRITIC_ITERS batches of data to train the
    critic, one batch to train the generator, and then one batch to evaluate
    one model performance.
    '''
    dataset = dataset.shuffle(num_data_points)
    dataset = dataset.batch(BATCH_SIZE)
    losses = []
    checkpoint_dir = save_dir + '/training_checkpoints'
    for epoch in range(epochs):
        batches = []
        count = 0
        for batch in dataset:
            count += 1
            batches.append(batch)
            if count % (NUM_CRITIC_ITERS + 2) == 0:
                wass_est = model.train_step(batches)
                print("Loss: {}".format(wass_est))
                losses.append(wass_est)
                batches.clear()
            if count % 700 == 0:
                print("Completed {} batches".format(count/7))


        if (epoch + 1) % 1 == 0:
            gen_filename = os.path.join(checkpoint_dir, 'gen_' + str(epoch + 1))
            critic_filename = os.path.join(checkpoint_dir, 'critic_' +
                    str(epoch + 1))
            print(gen_filename)
            model.generator.save_weights(gen_filename)
            model.critic.save_weights(critic_filename)

    return losses


def main():
    today = str(date.today())
    run_number = 0
    save_dir = '../../models/cWGAN/Run_' + today + '_' + str(run_number) 
    
    while (os.path.exists(save_dir)):
        run_number += 1
        save_dir = '../../models/cWGAN/Run_' + today + '_' + str(run_number) 

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert(os.path.isdir(save_dir))
    
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

    cwgan = model.cWGAN(NUM_CRITIC_ITERS, BATCH_SIZE)
    fname = os.path.join(save_dir, 'modelsummary.txt')
    with open(fname, 'w') as f:
        with redirect_stdout(f):
            cwgan.print_network()

    start = time.time()
    losses = train(cwgan, dataset, 1, num_data_points, save_dir)
    print("Time for {} epochs: {}".format(NUM_EPOCHS, time.time() - start))
    filename = save_dir + '/losses.txt'
    loss_df = pd.DataFrame({"Wass Est": pd.Series(losses)})
    loss_df.to_csv(filename, header='None', index='None', sep=' ')


if __name__ == "__main__":
    main()

