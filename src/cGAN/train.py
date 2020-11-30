import tensorflow as tf
from tensorflow import keras
import pandas as pd
import numpy as np
import model
from datetime import date
import time
from contextlib import redirect_stdout
import os


BATCH_SIZE=64
cross_entropy = tf.keras.losses.BinaryCrossentropy()
gen_optimizer = tf.keras.optimizers.Adam(3e-4)
disc_optimizer = tf.keras.optimizers.Adam(1e-4)
gen= model.make_generator()
disc= model.make_discriminator()


def discriminator_loss(real_output, fake_output):
    real_loss = cross_entropy(tf.ones_like(real_output), real_output)
    fake_loss = cross_entropy(tf.zeros_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    return total_loss

def generator_loss(fake_output):
    return cross_entropy(tf.ones_like(fake_output), fake_output)

@tf.function
def train_step(pJets, rJets):
    shape = tf.shape(pJets)
    noise = tf.random.normal(shape, 0, 1, tf.float32)

    with tf.GradientTape(persistent=True) as tape:
        generated_rJets = gen([pJets, noise], training=True)

        fake_output = disc([pJets, generated_rJets], training=True)
        real_output = disc([pJets, rJets], training=True)

        gen_loss = generator_loss(fake_output)
        disc_loss = discriminator_loss(real_output, fake_output)

    gen_grads = tape.gradient(gen_loss, gen.trainable_variables)
    disc_grads = tape.gradient(disc_loss, disc.trainable_variables)

    gen_optimizer.apply_gradients(zip(gen_grads, gen.trainable_variables))
    disc_optimizer.apply_gradients(zip(disc_grads,
        disc.trainable_variables))
    return gen_loss, disc_loss


@tf.function
def val_step(pJets, rJets):
    shape = tf.shape(pJets)
    noise = tf.random.normal(shape, 0, 1, tf.float32)

    generated_rJets = gen([pJets, noise], training=False)

    fake_output = disc([pJets, generated_rJets], training=False)
    real_output = disc([pJets, rJets], training=False)

    gen_loss = generator_loss(fake_output)
    disc_loss = discriminator_loss(real_output, fake_output)

    return gen_loss, disc_loss


#def train(train_dataset, val_pJets, val_rJets, epochs, save_dir):
def train(train_dataset, epochs, save_dir):
    checkpoint_dir = save_dir + '/training_checkpoints'
    train_gen_loss_list = []
    train_disc_loss_list= []

    #val_gen_loss_list= []
    #val_disc_loss_list= []

    for epoch in range(epochs):
        start = time.time()
        train_gen_loss_tot = 0
        train_disc_loss_tot = 0
        num_batch = 0

        for jet_batch in train_dataset:
            train_gen_loss, train_disc_loss = train_step(jet_batch[0], jet_batch[1])
            train_gen_loss_tot += train_gen_loss
            train_disc_loss_tot += train_disc_loss
            num_batch += 1

        #val_gen_loss, val_disc_loss = val_step(val_pJets, val_rJets)

        print('Time for epoch {} is {:.2f}. Gen Loss: {:.3f}  Disc Loss: {:.5f}'.format(epoch+1, time.time() - start, train_gen_loss_tot/num_batch, train_disc_loss_tot/num_batch))
        #print('Val Gen Loss: {:.3f}  Val Disc Loss: {:.3f}'.format(val_gen_loss, val_disc_loss))
        print("")

        train_gen_loss_list.append(train_gen_loss_tot.numpy()/num_batch)
        train_disc_loss_list.append(train_disc_loss_tot.numpy()/num_batch)
        #val_gen_loss_list.append(val_gen_loss.numpy())
        #val_disc_loss_list.append(val_disc_loss.numpy())

        if (epoch + 1) % 5 == 0:
            gen_filename = os.path.join(checkpoint_dir, 'gen_' + str(epoch + 1))
            disc_filename = os.path.join(checkpoint_dir, 'disc_' + str(epoch + 1))
            print(gen_filename)
            gen.save_weights(gen_filename)
            disc.save_weights(disc_filename)

    return (train_gen_loss_list, train_disc_loss_list)

    #return (train_gen_loss_list, train_disc_loss_list, val_gen_loss_list,
    #        val_disc_loss_list)


def main():
    today = str(date.today())
    run_number = 0
    save_dir = '../../models/cGAN/Run_' + today + '_' + str(run_number) 
    
    while (os.path.exists(save_dir)):
        run_number += 1
        save_dir = '../../models/cGAN/Run_' + today + '_' + str(run_number) 

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert(os.path.isdir(save_dir))
    
    
    fname = os.path.join(save_dir, 'modelsummary.txt')
    with open(fname, 'w') as f:
        with redirect_stdout(f):
            gen.summary()
            disc.summary()

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

    #index = int(0.8*len(data))
    #print("Number of training examples: {}".format(index))
    #print("Number of validation examples: {}".format(len(data) - index))
    #trainData = data[:index, :]
    #np.random.shuffle(trainData)
    np.random.shuffle(data)
    trainParton = data[:, :4]
    trainPf = data[:, 4:]
    train_dataset = tf.data.Dataset.from_tensor_slices((trainParton,
            trainPf))
    train_dataset = train_dataset.batch(BATCH_SIZE)
    #validate = data[index:, :]
    #validateParton = validate[:, :4]
    #validatePf = validate[:, 4:]

    
    #losses = train(train_dataset, validateParton, validatePf, 1000, save_dir)
    losses = train(train_dataset, 300, save_dir)
    loss_df = pd.DataFrame({'Gen Training Loss': pd.Series(losses[0]),
        'Disc Training Loss': pd.Series(losses[1])})#,
        #'Gen Validation Loss': pd.Series(losses[2]),
        #'Disc Validation Loss': pd.Series(losses[3])})

    filename = save_dir + '/losses.txt'
    loss_df.to_csv(filename, header='None', index='None', sep=' ')

if __name__ == "__main__":
    main()

