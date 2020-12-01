import tensorflow as tf
from tensorflow import keras
import pandas as pd
import numpy as np
import model
from datetime import date
import time
from contextlib import redirect_stdout
import os

def train(model, iterations, save_dir):
    losses = []
    checkpoint_dir = save_dir + '/training_checkpoints'
    for it in range(iterations):
        wass_est = model.train_step()
        print("Loss: {}".format(wass_est))
        losses.append(wass_est)

        if (it + 1) % 100 == 0:
            gen_filename = os.path.join(checkpoint_dir, 'gen_' + str(it + 1))
            critic_filename = os.path.join(checkpoint_dir, 'critic_' + str(it +
                1))
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
    num_batches = int(len(data)/64)
    print(num_batches)
    np.random.shuffle(data)
    trainParton = data[:, :4]
    trainPf = data[:, 4:]
    train_dataset = tf.data.Dataset.from_tensor_slices((trainParton,
            trainPf))

    cwgan = model.cWGAN(train_dataset, num_batches)
    fname = os.path.join(save_dir, 'modelsummary.txt')
    with open(fname, 'w') as f:
        with redirect_stdout(f):
            cwgan.print_network()

    losses = train(cwgan, 5000, save_dir)
    filename = save_dir + '/losses.txt'
    loss_df = pd.DataFrame({"Wass Est": pd.Series(losses)})
    loss_df.to_csv(filename, header='None', index='None', sep=' ')


if __name__ == "__main__":
    main()

