"""
Classes to implement and train a conditional wasserstein GAN on the parton data
as well as MNIST data. The cWGAN learns to generate samples y (either reco jet 
4-momenta or images of handwritten digits), conditioned onsome input x (either 
parton jet 4-momenta or a digit 0-9). 
"""


import tensorflow as tf
from tensorflow import keras
import numpy as np
import data_utils
import file_utils
import time
import os


class cWGAN:
    """
    Class implementing a conditional wasserstein generative adversarial
    network
    """

    def __init__(self, clip_value, noise_dims):
        """
        Constructor
        clip_value - value for weight clipping the critic
        noise_dims - dimension of the noise input to the generator
        """
        # hyper parameters recommended by paper
        self.clip_value = clip_value
        self.critic_optimizer = tf.keras.optimizers.RMSprop(lr=5e-5)
        self.generator_optimizer = tf.keras.optimizers.RMSprop(lr=5e-5)

        self.noise_dims = noise_dims
        self.generator = self.build_generator()
        self.critic = self.build_critic()

    def build_generator(self):
        noise = keras.Input(shape=(self.noise_dims,), name="noiseIn")
        x_in = keras.Input(shape=(4,), name="pjetIn")

        concat = keras.layers.concatenate([x_in, noise], name="concat")
        z = keras.layers.Dense(512, activation="relu")(concat)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        out = keras.layers.Dense(4, activation="relu")(z)

        return keras.Model([x_in, noise], out)

    def build_critic(self):
        x_in = keras.Input(shape=(4,))
        y_in = keras.Input(shape=(4,))
        concat = keras.layers.concatenate([x_in, y_in])



        z = keras.layers.Dense(512, activation="relu")(concat)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)

        out = keras.layers.Dense(1)(z)
        return keras.Model([x_in, y_in], out)

    @tf.function
    def critic_loss(self, real_output, fake_output):
        """
        The negative of the estimate of the wasserstein distance (negative
        because we want to perform gradient ascent on the critic)
        real_output - output of discriminator when given inputs x matched
        with real data y jets
        fake_output - output of discriminator when given inputs x and outputs y
        from the generator
        returns - wasserstein loss
        """
        loss = -(tf.math.reduce_mean(real_output) - tf.math.reduce_mean(fake_output))
        return loss

    @tf.function
    def generator_loss(self, fake_output):
        """
        Estimate of wasserstein loss for the generator
        fake_output - output of discriminator when given inputs x and outputs y
        from the generator
        returns - wasserstein loss
        """

        loss = -tf.math.reduce_mean(fake_output)
        return loss

    def clip_critic_weights(self):
        """
        Clip the weights of the critic to value set by self.clip_value
        """
        for l in self.critic.layers:
            new_weights = []
            for i in range(len(l.weights)):
                new_weights.append(
                    tf.clip_by_value(l.weights[i], -self.clip_value, self.clip_value)
                )
            l.set_weights(new_weights)

    @tf.function
    def train_critic(self, x, y):
        """
        Train critic on one batch of data
        x - batch of input data
        y - batch of matching output data
        returns - the critic loss for the batch
        """
        noise = tf.random.uniform(
            (tf.shape(x)[0], self.noise_dims), 0, 1, tf.float32
        )
        with tf.GradientTape(persistent=True) as tape:
            predicted_y = self.generator([x, noise], training=False)
            real_output = self.critic([x, y], training=True)
            fake_output = self.critic([x, predicted_y], training=True)

            critic_loss_val = self.critic_loss(real_output, fake_output)

        critic_grads = tape.gradient(critic_loss_val, self.critic.trainable_variables)

        self.critic_optimizer.apply_gradients(
            zip(critic_grads, self.critic.trainable_variables)
        )

        return critic_loss_val

    @tf.function
    def train_generator(self, x):
        """
        Train generator on one batch of data
        x - batch of input data
        """
        noise = tf.random.uniform(
            (tf.shape(x)[0], self.noise_dims), 0, 1, tf.float32
        )

        with tf.GradientTape() as tape:
            generated_rJets = self.generator([x, noise], training=True)
            fake_output = self.critic([x, generated_rJets], training=False)
            generator_loss_val = self.generator_loss(fake_output)

        generator_grads = tape.gradient(
            generator_loss_val, self.generator.trainable_variables
        )
        self.generator_optimizer.apply_gradients(
            zip(generator_grads, self.generator.trainable_variables)
        )

        return generator_loss_val

    @tf.function
    def make_generator_predictions(self, x):
        noise = tf.random.uniform(
            (tf.shape(x)[0], self.noise_dims), 0, 1, tf.float32
        )
        predictions = self.generator([x, noise], training=False)
        return predictions


class Trainer:
    def __init__(self, params_dict):
        self.save_dir = file_utils.make_save_directory("cWGAN")
        self.num_critic_iters = params_dict["num_critic_iters"]
        self.batch_size = params_dict["batch_size"]
        clip_value = params_dict["clip_value"]
        noise_dims = params_dict["noise_dims"]
        self.model = cWGAN(clip_value, noise_dims)

        self.data = data_utils.load_jet_data(params_dict["data_path"])
        self.epochs = params_dict["epochs"]
        self.num_training_examples = len(self.data[0])
        self.weight_saving_interval = params_dict["weight_saving_interval"]

        self.critic_losses = []
        self.generator_losses = []
        self.wass_estimates = []

    def sample_batch_of_data(self):
        indices = np.random.choice(
            np.arange(self.num_training_examples), self.batch_size, replace=False
        )
        labels = self.data[0][indices]
        images = self.data[1][indices]
        return labels, images

    def take_critic_step(self):
        x, y = self.sample_batch_of_data()
        critic_loss = self.model.train_critic(x, y)
        self.critic_losses.append(critic_loss)
        self.model.clip_critic_weights()

    def take_generator_step(self):
        x, y = self.sample_batch_of_data()
        generator_loss = self.model.train_generator(x)
        self.generator_losses.append(generator_loss)

        predicted_y = self.model.make_generator_predictions(x)
        real_output = self.model.critic([x, y], training=False)
        fake_output = self.model.critic([x, predicted_y], training=False)

        wass_estimate = -self.model.critic_loss(real_output, fake_output)
        self.wass_estimates.append(wass_estimate)


    def train(self):
        """
        Training loop for the cWGAN. An "epoch" is considered to be when the generator
        has seen the same number of examples as are in the data set (note that since the
        batches are randomly sampled, it won't actually get trained on all the data
        each epoch).
        """
        batches_per_epoch = self.num_training_examples // self.batch_size
        for epoch in range(self.epochs):
            start = time.time()
            for batch_number in range(5):#range(batches_per_epoch):
                # train critic for num_critic_iters
                for critic_iter in range(self.num_critic_iters):
                    self.take_critic_step()
                # train generator
                self.take_generator_step()
                   
                iteration = epoch * batches_per_epoch + batch_number
                if iteration % self.weight_saving_interval == 0:
                    self.save_weights(iteration)
                print(
                    "Iteration: {}  Wasserstein Estimate: {}".format(
                        iteration, self.wass_estimates[-1]
                    )
                )
            print("Time for epoch {}: {:1f}s".format(epoch, time.time() - start))

    def save_weights(self, iteration):
        checkpoint_dir = self.save_dir + "/training_checkpoints"
        gen_filename = os.path.join(checkpoint_dir, "gen_" + str(iteration))
        critic_filename = os.path.join(checkpoint_dir, "critic_" + str(iteration))
        print("Saving generator weights at {}".format(gen_filename))
        self.model.generator.save_weights(gen_filename)
        self.model.critic.save_weights(critic_filename)

    def save_losses(self):
        critic_loss_dict = {
            "Critic Loss": self.critic_losses,
        }
        file_utils.save_losses(self.save_dir, critic_loss_dict, "critic_")
        generator_loss_dict = {
            "Generator Loss": self.generator_losses,
            "Wasserstein Estimates": self.wass_estimates,
        }
        file_utils.save_losses(self.save_dir, generator_loss_dict, "generator_")

    def save_model(self):
        file_utils.save_network(self.save_dir, model_path="./cWGAN.py")

    def save_params(self, params_dict):
        file_utils.save_params(self.save_dir, params_dict)

class cWGAN_mnist(cWGAN):
    """
    Subclass of cWGAN class to be used with MNIST data.
    """

    def build_generator(self):
        """
        Override the build generator method from parent class. Instead returns
        CNN based generator to use with MNIST data.
        """
        noise = keras.Input(shape=(self.noise_dims,))
        number_input = keras.Input(shape=(10,))

        input1 = keras.layers.Dense(10, activation="relu")(number_input)
        input1 = keras.layers.Dense(32, activation="relu")(input1)

        input2 = keras.layers.Dense(self.noise_dims, activation="relu")(noise)
        input2 = keras.layers.Dense(self.noise_dims, activation="relu")(input2)

        concat = keras.layers.concatenate([input1, input2])
        out = keras.layers.Dense(7 * 7 * 256, activation="relu")(concat)

        out = keras.layers.Reshape((7, 7, 256))(out)
        out = keras.layers.Conv2DTranspose(
            256, (5, 5), strides=(1, 1), padding="same", use_bias=False
        )(out)
        out = keras.layers.LeakyReLU()(out)
        out = keras.layers.Conv2DTranspose(
            128, (5, 5), strides=(2, 2), padding="same", use_bias=False
        )(out)
        out = keras.layers.LeakyReLU()(out)
        out = keras.layers.Conv2DTranspose(
            1, (5, 5), strides=(2, 2), padding="same", use_bias=False, activation="tanh"
        )(out)

        return keras.Model([number_input, noise], out)

    def build_critic(self):
        """
        Override the build discriminator method from parent class. Instead returns
        CNN based discriminator to use with MNIST data.
        """

        image = keras.Input(shape=(28, 28, 11))

        z = keras.layers.Conv2D(64, (5, 5), strides=(2, 2), padding="same")(image)
        z = keras.layers.LeakyReLU()(z)
        z = keras.layers.Conv2D(128, (5, 5), strides=(2, 2), padding="same")(z)
        z = keras.layers.LeakyReLU()(z)
        z = keras.layers.Conv2D(256, (5, 5), strides=(2, 2), padding="same")(z)
        z = keras.layers.LeakyReLU()(z)
        z = keras.layers.Conv2D(512, (5, 5), strides=(2, 2), padding="same")(z)
        z = keras.layers.LeakyReLU()(z)
        out = keras.layers.Conv2D(1, 2, 1)(z)

        return keras.Model(image, out)

    @tf.function
    def train_critic(self, labels, images):
        """
        Train critic on one batch of data
        labels - batch of parton jet 4-momenta
        images - batch of maching reco jet 4-momenta
        returns - the critic loss for the batch
        """
        noise = tf.random.uniform(
            (tf.shape(labels)[0], self.noise_dims), 0, 1, tf.float32
        )
        concat_real = data_utils.concatenate_images_labels(images, labels)
        with tf.GradientTape(persistent=True) as tape:
            generated_images = self.generator([labels, noise], training=False)
            concat_fake = data_utils.concatenate_images_labels(generated_images, labels)
            real_output = self.critic(concat_real, training=True)
            fake_output = self.critic(concat_fake, training=True)

            critic_loss_val = self.critic_loss(real_output, fake_output)

        critic_grads = tape.gradient(critic_loss_val, self.critic.trainable_variables)

        self.critic_optimizer.apply_gradients(
            zip(critic_grads, self.critic.trainable_variables)
        )

        return critic_loss_val

    @tf.function
    def train_generator(self, labels):
        """
        Train generator on one batch of data
        """
        noise = tf.random.uniform(
            (tf.shape(labels)[0], self.noise_dims), 0, 1, tf.float32
        )

        with tf.GradientTape() as tape:
            generated_images = self.generator([labels, noise], training=True)
            concat_fake = data_utils.concatenate_images_labels(generated_images, labels)
            fake_output = self.critic(concat_fake, training=False)
            generator_loss_val = self.generator_loss(fake_output)

        generator_grads = tape.gradient(
            generator_loss_val, self.generator.trainable_variables
        )
        self.generator_optimizer.apply_gradients(
            zip(generator_grads, self.generator.trainable_variables)
        )

        return generator_loss_val


class MNISTTrainer(Trainer):
    def __init__(self, params_dict):
        self.save_dir = file_utils.make_save_directory("cWGAN")
        self.num_critic_iters = params_dict["num_critic_iters"]
        self.batch_size = params_dict["batch_size"]
        clip_value = params_dict["clip_value"]
        noise_dims = params_dict["noise_dims"]
        self.model = cWGAN_mnist(clip_value, noise_dims)

        self.data = data_utils.load_mnist_data()
        self.epochs = params_dict["epochs"]
        self.num_training_examples = len(self.data[0])
        self.weight_saving_interval = params_dict["weight_saving_interval"]

        self.critic_losses = []
        self.generator_losses = []
        self.wass_estimates = []

    
    def take_generator_step(self):
        """
        Override function from parent class to handle the image - label
        concatenation necessary for the critic
        """
        labels, images = self.sample_batch_of_data()
        concat_real = data_utils.concatenate_images_labels(images, labels)
        generator_loss = self.model.train_generator(labels)
        predicted_images = self.model.make_generator_predictions(labels)
        concat_fake = data_utils.concatenate_images_labels(
            predicted_images, labels
        )
        self.generator_losses.append(generator_loss)
        real_output = self.model.critic(concat_real, training=False)
        fake_output = self.model.critic(concat_fake, training=False)
        wass_estimate = -self.model.critic_loss(real_output, fake_output)
        self.wass_estimates.append(wass_estimate)

    """
    def train(self):
        
        Training loop for the cWGAN. An "epoch" is considered to be when the generator
        has seen the same number of examples as are in the data set (note that since the
        batches are randomly sampled, it won't actually get trained on all the data
        each epoch).
        
        batches_per_epoch = self.num_training_examples // self.batch_size
        for epoch in range(self.epochs):
            start = time.time()
            for batch_number in range(batches_per_epoch):
                # train critic for num_critic_iters
                for critic_iter in range(self.num_critic_iters):
                    labels, images = self.sample_batch_of_data()
                    critic_loss = self.model.train_critic(labels, images)
                    self.critic_losses.append(critic_loss)
                    self.model.clip_critic_weights()
                # train generator
                
                iteration = epoch * batches_per_epoch + batch_number
                if iteration % self.weight_saving_interval == 0:
                    self.save_weights(iteration)
                print(
                    "Iteration: {}  Wasserstein Estimate: {}".format(
                        iteration, wass_estimate
                    )
                )
            print("Time for epoch {}: {:1f}s".format(epoch, time.time() - start))
    """


def main():

    net = cWGAN(0.1, 100)
    net.generator.summary()
    net.critic.summary()


if __name__ == "__main__":
    main()
