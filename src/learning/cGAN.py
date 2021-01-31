import os
import tensorflow as tf
from tensorflow import keras
import file_utils
import data_utils
import time
import numpy as np


class cGAN:
    """Class implementing a conditional generative adversarial network described
    here: https://arxiv.org/pdf/1411.1784.pdf. The generator G takes in noise z,
    some conditional data x, and produces an output y. The network tries to learn
    to produce samples from the conditional distribution P(y|x).
    """

    def __init__(self, noise_dims, gen_lr, disc_lr):
        """Constructor

        Args:
            noise_dims (int): Size of noise vector input to generator
            gen_lr (float): Learning rate for generator optimizer
            disc_lr (float): Learning rate for discriminator optimizer
        """
        self.noise_dims = noise_dims

        self.generator = self.make_generator()
        self.discriminator = self.make_discriminator()

        self.gen_optimizer = tf.keras.optimizers.Adam(gen_lr)
        self.discriminator_optimizer = tf.keras.optimizers.Adam(disc_lr)

        self.cross_entropy = tf.keras.losses.BinaryCrossentropy()

    def make_generator(self):
        noise = keras.Input(shape=(self.noise_dims,), name="noiseIn")
        pJet = keras.Input(shape=(4,), name="pjetIn")

        x = keras.layers.Dense(32, activation="relu")(pJet)
        x = keras.layers.Dense(32, activation="relu")(x)
        x = keras.layers.Dense(32, activation="relu")(x)

        y = keras.layers.Dense(32, activation="relu")(noise)
        y = keras.layers.Dense(32, activation="relu")(y)
        y = keras.layers.Dense(32, activation="relu")(y)

        concat = keras.layers.concatenate([x, y])
        out = keras.layers.Dense(32, activation="relu")(concat)
        out = keras.layers.Dense(32, activation="relu")(concat)
        out = keras.layers.Dense(4, name="out")(out)
        return keras.Model([pJet, noise], out)

    def make_discriminator(self):
        pJet = keras.Input(shape=(4,))
        rJet = keras.Input(shape=(4,))

        x = keras.layers.Dense(32, activation="relu")(pJet)
        x = keras.layers.Dense(32, activation="relu")(x)
        x = keras.layers.Dense(32, activation="relu")(x)

        y = keras.layers.Dense(32, activation="relu")(rJet)
        y = keras.layers.Dense(32, activation="relu")(y)
        y = keras.layers.Dense(32, activation="relu")(y)

        concat = keras.layers.concatenate([x, y])
        out = keras.layers.Dense(32, activation="relu")(concat)
        out = keras.layers.Dense(32, activation="relu")(concat)
        out = keras.layers.Dense(1, activation="sigmoid")(out)
        return keras.Model([pJet, rJet], out)

    def discriminator_loss(self, real_output, fake_output):
        """Calculate binary crossentropy loss for the discriminator

        Args:
            real_output (tf.Tensor): Output from discriminator after being given parton jets matched with real reco jets
            fake_output (tf.Tensor): Output from discriminator after being given parton jets matched with reco jets from the

        Returns:
            tf.Tensor: loss for the discriminator
        """

        real_loss = self.cross_entropy(tf.ones_like(real_output), real_output)
        fake_loss = self.cross_entropy(tf.zeros_like(fake_output), fake_output)
        total_loss = real_loss + fake_loss
        return total_loss

    def generator_loss(self, fake_output):
        """Calculate binary crossentropy loss for the generator

        Args:
            fake_output (tf.Tensor): Output from the discriminator after being given parton jets and reco jets from the generator

        Returns:
            tf.Tensor: Loss for the generator
        """

        return self.cross_entropy(tf.ones_like(fake_output), fake_output)

    def train_step(self, x, y):
        """One forward pass and backpropagation pass for a batch of data

        Args:
            x (tf.Tensor, shape(batch_size, 4)): Data the network is conditioned on
            y (tf.Tensor, shape(batch_size, 4)): [description]

        Returns:
            gen_loss (tf.Tensor): The generator loss for the batch
            disc_loss (tf.Tensor): The discriminator loss for the batch
        """
        noise = tf.random.uniform((tf.shape(x)[0], self.noise_dims), 0, 1, tf.float32)

        with tf.GradientTape(persistent=True) as tape:
            generated_y = self.generator([x, noise], training=True)

            fake_output = self.discriminator([x, generated_y], training=True)
            real_output = self.discriminator([x, y], training=True)

            gen_loss = self.generator_loss(fake_output)
            disc_loss = self.discriminator_loss(real_output, fake_output)

        gen_grads = tape.gradient(gen_loss, self.generator.trainable_variables)
        disc_grads = tape.gradient(disc_loss, self.discriminator.trainable_variables)

        self.gen_optimizer.apply_gradients(
            zip(gen_grads, self.generator.trainable_variables)
        )
        self.discriminator_optimizer.apply_gradients(
            zip(disc_grads, self.discriminator.trainable_variables)
        )

        return gen_loss, disc_loss


class Trainer:
    def __init__(self, params_dict):
        """Constructor

        Args:
            params_dict (dict): Dictionary with model hyperparams
        """
        self.save_dir = file_utils.make_save_directory("cGAN")
        self.batch_size = params_dict["batch_size"]
        noise_dims = params_dict["noise_dims"]
        gen_lr = params_dict["gen_lr"]
        disc_lr = params_dict["disc_lr"]
        self.model = cGAN(noise_dims, gen_lr, disc_lr)

        dataset = data_utils.load_jet_data(params_dict["data_path"])
        tf_dataset = tf.data.Dataset.from_tensor_slices((dataset[0], dataset[1]))
        self.data = tf_dataset.batch(self.batch_size)
        self.epochs = params_dict["epochs"]

        self.discriminator_losses = []
        self.generator_losses = []
        self.weight_saving_interval = params_dict["weight_saving_interval"]

    def train(self):
        for epoch in range(self.epochs):
            start = time.time()

            total_generator_loss = 0
            total_discriminator_loss = 0

            num_batches = 0

            for batch in self.data:
                generator_loss, discriminator_loss = self.model.train_step(
                    batch[0], batch[1]
                )
                total_generator_loss += generator_loss
                total_discriminator_loss += discriminator_loss
                num_batches += 1

            avg_generator_loss = total_generator_loss.numpy() / num_batches
            avg_discriminator_loss = total_discriminator_loss.numpy() / num_batches

            self.discriminator_losses.append(avg_discriminator_loss)
            self.generator_losses.append(avg_generator_loss)

            print(
                "Time for epoch {} is {:.2f}. Gen Loss: {:.3f}   Disc Loss {:.3f}".format(
                    epoch,
                    time.time() - start,
                    avg_generator_loss,
                    avg_discriminator_loss,
                )
            )

            if epoch % self.weight_saving_interval == 0:
                self.save_weights(epoch)

    def save_weights(self, epoch):
        """Save the weights of the model to the save directory.

        Args:
            epoch (int): Training epoch
        """
        checkpoint_dir = self.save_dir + "/training_checkpoints"
        gen_filename = os.path.join(checkpoint_dir, "gen_" + str(epoch))
        print("Saving generator weights at {}".format(gen_filename))
        self.model.generator.save_weights(gen_filename)

    def save_losses(self):
        """Save losses to txt file
        """
        loss_dict = {
            "Gen Loss": self.generator_losses,
            "Discriminator Loss": self.discriminator_losses,
        }
        file_utils.save_losses(self.save_dir, loss_dict)

    def save_model(self):
        """Copy cGAN.py to save dir, for later model evaluation
        """
        file_utils.save_network(self.save_dir, model_path="./cGAN.py")

    def save_params(self, params_dict):
        """Save model parameters

        Args:
            params_dict (dict): Contains training/model parameters
        """
        file_utils.save_params(self.save_dir, params_dict)


def main():
    pass


if __name__ == "__main__":
    main()
