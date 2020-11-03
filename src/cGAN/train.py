import tensorflow as tf
from tensorflow import keras
import numpy as np
import model
from datetime import date
import time
from contextlib import redirect_stdout
import os

BATCH_SIZE=64


def main():
    '''
    today = str(date.today())
    run_number = 0
    save_dir = '../models/cGAN/Run_' + str(run_number) + '_' + today
    
    while (os.path.exists(save_dir)):
        run_number += 1
        save_dir = '../models/cGAN/Run_' + str(run_number)  + '_' + today

    print("SAVE DIR: " + save_dir)
    os.makedirs(save_dir)
    assert(os.path.isdir(save_dir))
    '''
    
    gen= model.make_generator()
    disc= model.make_discriminator()

    gen_optimizer = tf.keras.optimizers.Adam(1e-4)
    disc_optimizer = tf.keras.optimizers.Adam(1e-4)
    '''
    fname = os.path.join(save_dir, 'modelsummary.txt')
    with open(fname, 'w') as f:
        with redirect_stdout(f):
            gen.summary()
            disc.summary()

    '''
    data = np.loadtxt('../../data/txt/matchttbarDataTot.txt', skiprows=2)
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

    index = int(0.8*len(data))
    print("Number of training examples: {}".format(index))
    print("Number of validation examples: {}".format(len(data) - index))
    train = data[:index, :]
    np.random.shuffle(train)
    trainParton = train[:, :4]
    trainPf = train[:, 4:]
    train_dataset = tf.data.Dataset.from_tensor_slices((trainParton,
            trainPf))
    train_dataset.batch(BATCH_SIZE)
    validate = data[index:, :]
    validateParton = validate[:, :4]
    validatePf = validate[:, 4:]



    cross_entropy = tf.keras.losses.BinaryCrossentropy()

    def discriminator_loss(real_output, fake_output):
        real_loss = cross_entropy(tf.ones_like(real_output), real_output)
        fake_loss = cross_entropy(tf.zeros_like(fake_output), fake_output)
        total_loss = real_loss + fake_loss
        return total_loss

    def generator_loss(fake_output):
        return cross_entropy(tf.ones_like(fake_output), fake_output)

    def train_step(pJets, rJets):
        noise = tf.random.normal([BATCH_SIZE, 4], 0, 1, tf.float32)

        with tf.GradientTape() as gen_tape, tf.GradientTape() as disc_tape:
            generated_rJets = gen([pJets, noise], training=True)

            real_output = disc([pJets, generated_rJets], training=True)
            fake_output = disc([pJets, rJets], training=True)

            gen_loss = generator_loss(fake_output)
            disc_loss = discriminator_loss(real_output, fake_output)

        gen_grads = gen_tape.gradient(gen_loss, gen.trainable_variables)
        disc_grads = disc_tape.gradient(disc_loss, disc.trainable_variables)

        gen_optimizer.apply_gradients(zip(gen_grads,
            gen.trainable_variables))
        disc_optimizer.apply_gradiets(zip(disc_grads,
            disc.trainable_variables))

    def train(dataset, epochs):
        for epoch in range(epochs):
            start = time.time()

            for jet_batch in dataset:
                train_step(jet_batch[0], jet_batch[1])

        print('Time for epoch {} is {}'.format(epoch+1, time.time() - start))

    train(train_dataset, 10)

            


if __name__ == "__main__":
    main()

