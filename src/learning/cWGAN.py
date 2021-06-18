"""
Classes to implement and train a conditional wasserstein GAN on the parton data
as well as MNIST data. The cWGAN learns to generate samples y (either reco jet 
4-momenta or images of handwritten digits), conditioned on some input x (either 
parton jet 4-momenta or a digit 0-9). 
"""


import tensorflow as tf
from tensorflow import keras
import numpy as np
import data_utils
import file_utils
import time
import os
import sys


class cWGAN:
    """Class implementing a conditional wasserstein generative adversarial network"""

    def __init__(self, noise_dims, optimizer, gen_lr, critic_lr, gp_weight, load_previous, weights_path, iteration):
        """Constructor

        Args:
            noise_dims (int): Size of noise vector input to the generator
            optimizer (str): Choice of optimizer - either RMSprop or Adam
            gen_lr (float): learning rate of the generator optimizer
            critic_lr (float): learning rate of the critic optimizer
            gp_weight (float): Weight for the gradient penalty in the loss function
            load_previous (bool): Whether to load weights from a previous run
            weights_path (path-like): If load_previous is true, path to model weights
            iteration (int): If load_previous is true, the iteration of weights to load
        """
        if optimizer == "RMSprop":
            self.critic_optimizer = tf.keras.optimizers.RMSprop(lr=critic_lr)
            self.generator_optimizer = tf.keras.optimizers.RMSprop(lr=gen_lr)
        elif optimizer == "Adam":
            self.critic_optimizer = tf.keras.optimizers.Adam(lr=critic_lr)
            self.generator_optimizer = tf.keras.optimizers.Adam(lr=gen_lr)
        else:
            print("There was an error choosing an optimizer")
            print("Please specify one of RMSprop or Adam in your configuration file")
            sys.exit(1)

        self.gp_weight = gp_weight

        self.noise_dims = noise_dims
        self.generator = self.build_generator()
        self.critic = self.build_critic()

        if load_previous:
            self.initialize_model(weights_path, iteration)

    def initialize_model(self, weights_path, iteration):
        """Load weights from previous run"""
        self.generator.load_weights(weights_path + 'gen_' + str(iteration))
        self.critic.load_weights(weights_path + 'critic_' + str(iteration))


    def build_generator(self):
        noise = keras.Input(shape=(self.noise_dims,), name="noiseIn")
        x_in = keras.Input(shape=(12,), name="pjetIn")

        concat = keras.layers.concatenate([x_in, noise], name="concat")
        z = keras.layers.Dense(512, activation="relu")(concat)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        z = keras.layers.Dense(512, activation="relu")(z)
        out = keras.layers.Dense(12)(z)

        return keras.Model([x_in, noise], out)

    def build_critic(self):
        x_in = keras.Input(shape=(12,))
        y_in = keras.Input(shape=(12,))
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
        """Calculates the negative of the wasserstein distance (negative because we
        want to perform gradient ascent - not descent - on the critic) between the
        target distribution, and the distribution generated by the generator.

        Args:
            real_output (tf.Tensor): Output of the discriminator when given inputs x
            matched with real data y
            fake_output (tf.Tensor): Output of the discriminator when given inputs x
            and outputs y from the generator

        Returns:
            tf.Tensor: Estimate of Wasserstein distance
        """

        loss = -(tf.math.reduce_mean(real_output) - tf.math.reduce_mean(fake_output))
        return loss

    @tf.function
    def generator_loss(self, fake_output):
        """Estimate the Wasserstein loss for the generator

        Returns:
            tf.Tensor: Wasserstein loss
        """

        loss = -tf.math.reduce_mean(fake_output)
        return loss

    @tf.function
    def interpolate_data(self, x_real, x_gen, y_real, y_gen):
        """Interpolate between data points as described here: https://arxiv.org/pdf/1704.00028.pdf

        The gradient penalty acts on data sampled from straight lines between points in
        the real distribution P_r and the generated distribution P_g, so for points
        (x_real, y_real) ~ P_r and (x_gen, y_gen) ~ P_g, we want 

        (x_new, y_new) = t*(x_gen, y_gen) + (1-t)*(x_real, y_real)
              = t*((x_gen, y_gen) - (x_real, y_real)) + (x_real, y_real)

        where t is in (0, 1).
        Args:
            x_real (tf.Tensor): Batch of input data sampled from the real distribution
            x_gen (tf.Tensor): Batch of input data sampled from the generated distribution
            y_real (tf.Tensor): Batch of data sampled from real distribution
            y_gen (tf.Tensor): Batch of data sampled from fake distribution

        Returns:
            tf.Tensor: Interpolated batch of data
        """
        batch_size = tf.shape(y_real)[0]
        t = tf.random.normal([batch_size, 1], 0, 1, tf.float32)
        y_diff = y_gen - y_real
        x_diff = x_gen - x_real
        y_new = t * y_diff + y_real
        x_new = t * x_diff + x_real

        return x_new, y_new

    def gradient_penalty(self, x_real, x_gen, y_real, y_gen):
        """Calculate the gradient penalty. See here for explantion: https://arxiv.org/pdf/1704.00028.pdf

        Args:
            x_real (tf.Tensor): Batch of real data that the distribution is conditioned on
            x_gen (tf.Tensor): Batch of gen data that the distribution is conditioned on
            y_real (tf.Tensor): Batch of data sampled from real distribution
            y_gen (tf.Tensor): Batch of data sampled from fake distribution

        Returns:
            tf.Tensor: The gradient penalty
        """
        x_interpolated, y_interpolated = self.interpolate_data(
            x_real, x_gen, y_real, y_gen
        )

        with tf.GradientTape() as gp_tape:
            gp_tape.watch(y_interpolated)
            gp_tape.watch(x_interpolated)
            pred = self.critic([x_interpolated, y_interpolated], training=True)

        grads = gp_tape.gradient(pred, [x_interpolated, y_interpolated])
        concat_grads = tf.concat([grads[0], grads[1]], 0)
        norm = tf.sqrt(tf.reduce_sum(tf.square(concat_grads), axis=1))
        gp = tf.reduce_mean((norm - 1.0) ** 2)
        return gp

    @tf.function
    def train_critic(self, x1, y1, x2, y2):
        """Train critic on one batch of data.

        Note that training the critic requires sampling two batches of data from the 
        joint distribution in order to calculate the gradient penalty.

        Args:
            x1 (tf.Tensor): First batch of input data generator is conditioned on
            y1 (tf.Tensor): First batch of corresponding real output data
            x2 (tf.Tensor): Second batch of input data generator is conditioned on
            y2 (tf.Tensor): Second batch of corresponding real output data

        Returns:
            tf.Tensor: Critic loss for the batch
        """

        noise = tf.random.uniform((tf.shape(x1)[0], self.noise_dims), 0, 1, tf.float32)
        with tf.GradientTape(persistent=True) as tape:
            predicted_y = self.generator([x2, noise], training=False)
            real_output = self.critic([x2, y2], training=True)
            fake_output = self.critic([x2, predicted_y], training=True)

            critic_loss_val = self.critic_loss(
                real_output, fake_output
            ) + self.gp_weight * self.gradient_penalty(x1, x2, y1, predicted_y)

        critic_grads = tape.gradient(critic_loss_val, self.critic.trainable_variables)

        self.critic_optimizer.apply_gradients(
            zip(critic_grads, self.critic.trainable_variables)
        )

        return critic_loss_val

    @tf.function
    def train_generator(self, x):
        """Train generator on one batch of data

        Args:
            x (tf.Tensor): Batch of data generator is conditioned on

        Returns:
            tf.Tensor: The loss for the batch
        """

        noise = tf.random.uniform((tf.shape(x)[0], self.noise_dims), 0, 1, tf.float32)

        with tf.GradientTape() as tape:
            predicted_y = self.generator([x, noise], training=True)
            fake_output = self.critic([x, predicted_y], training=False)
            generator_loss = self.generator_loss(fake_output)

        generator_grads = tape.gradient(
            generator_loss, self.generator.trainable_variables
        )
        self.generator_optimizer.apply_gradients(
            zip(generator_grads, self.generator.trainable_variables)
        )

        return generator_loss

    @tf.function
    def make_generator_predictions(self, x):
        """Generate predictions from the generator without training.

        Args:
            x (tf.Tensor): Input data

        Returns:
            tf.Tensor: Generator output
        """
        noise = tf.random.uniform((tf.shape(x)[0], self.noise_dims), 0, 1, tf.float32)
        predictions = self.generator([x, noise], training=False)
        return predictions


class Trainer:
    """Class used to train the cWGAN"""

    def __init__(self, params_dict):
        """Constructor

        Args:
            params_dict (dict): Contains training/model parameters
        """
        self.save_dir = file_utils.make_save_directory("cWGAN")
        self.num_critic_iters = params_dict["num_critic_iters"]
        self.batch_size = params_dict["batch_size"]
        gen_lr = params_dict["gen_lr"]
        critic_lr = params_dict["critic_lr"]
        optimizer = params_dict["optimizer"]
        noise_dims = params_dict["noise_dims"]
        gp_weight = params_dict["gp_weight"]
        load_previous = params_dict["load_previous"]
        weights_path = params_dict["weights_path"]
        iteration = params_dict["iteration"]

        self.model = cWGAN(
            noise_dims, optimizer, gen_lr, critic_lr, gp_weight, load_previous, weights_path, iteration
        )

        
        data_path = params_dict["data_path"]
        if params_dict["data_type"] == "angular":
            if params_dict["data_scaling"] == "inverse":
                self.data = data_utils.load_jet_data_inverse_scaling(data_path)
            elif params_dict["data_scaling"] == "log":
                self.data = data_utils.load_jet_data_log_scaling(data_path)
            elif params_dict["data_scaling"] == "minmax":
                self.data = data_utils.load_jet_data(data_path)
            else:
                print("There was an error loading the data")
                print(
                    "Please specify a data scaling scheme in your configuration file: The options are: inverse, log, and minmax"
                )
                sys.exit(1)
        elif params_dict["data_type"] == "cartesian":
            self.data = data_utils.load_jet_data_log_scaling_cartesian(data_path)
        else:
            print("There was an error loading the data")
            print("Please specify the form of the data your are training on, either:")
            print("    angular: Four-momenta in form (Pt, Eta, Phi, E)")
            print("    cartesian: Four-momenta in form (Px, Py, Pz, E)")
            sys.exit(1)

        self.epochs = params_dict["epochs"]
        self.num_training_examples = len(self.data[0])
        self.weight_saving_interval = params_dict["weight_saving_interval"]

        self.critic_losses = []
        self.generator_losses = []
        self.wass_estimates = []



    def sample_batch_of_data(self):
        """Randomly sample self.batch_size (x, y) pairs from the dataset

        Returns:
            tf.Tensor: Batch of data
        """
        indices = np.random.choice(
            np.arange(self.num_training_examples), self.batch_size, replace=False
        )
        x = self.data[0][indices]
        y = self.data[1][indices]
        return x, y

    def take_critic_step(self):
        """Sample a batch of data and do one forward pass and backpropagation step
        for the critic
        """
        x1, y1 = self.sample_batch_of_data()
        x2, y2 = self.sample_batch_of_data()
        critic_loss = self.model.train_critic(x1, y1, x2, y2)
        self.critic_losses.append(critic_loss)

    def take_generator_step(self):
        """Sample a batch of data and do one forward pass and backpropagation step
        for the generator. Make an estimate of the wasserstein distance between the
        generator post backprop step and the target distribution.
        """
        x, y = self.sample_batch_of_data()
        generator_loss = self.model.train_generator(x)
        self.generator_losses.append(generator_loss)

        predicted_y = self.model.make_generator_predictions(x)
        real_output = self.model.critic([x, y], training=False)
        fake_output = self.model.critic([x, predicted_y], training=False)

        wass_estimate = -self.model.critic_loss(real_output, fake_output)
        self.wass_estimates.append(wass_estimate)

    def train(self):
        """Training loop for the cWGAN.

        An "epoch" is considered to be when the generator has seen the same number of
        examples as are in the data set (note that since the batches are randomly
        sampled, it won't actually get trained on all the data each epoch).
        """
        batches_per_epoch = self.num_training_examples // self.batch_size
        for epoch in range(self.epochs):
            start = time.time()
            for batch_number in range(batches_per_epoch):
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
        """Save weights of the model to the save directory

        Args:
            iteration (int): Current generator iteration
        """
        checkpoint_dir = self.save_dir + "/training_checkpoints"
        gen_filename = os.path.join(checkpoint_dir, "gen_" + str(iteration))
        critic_filename = os.path.join(checkpoint_dir, "critic_" + str(iteration))
        print("Saving generator weights at {}".format(gen_filename))
        print("Saving generator weights at {}".format(critic_filename))
        self.model.generator.save_weights(gen_filename)
        self.model.critic.save_weights(critic_filename)

    def save_losses(self):
        """Save training losses to a txt file"""
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
        """Copy cWGAN.py to save dir, for later model evaluation"""
        file_utils.save_network(self.save_dir, model_path="./cWGAN.py")

    def save_params(self, params_dict):
        """Save model parameters

        Args:
            params_dict (dict): Contains training/model parameters
        """
        file_utils.save_params(self.save_dir, params_dict)


class cWGAN_mnist(cWGAN):
    """Subclass of cWGAN class to be used with MNIST data."""

    def build_generator(self):
        """Override the build generator method from parent class. Instead returns
        CNN based generator to use with MNIST data.

        Returns:
            keras.Model: Generator model
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
        """Override the build discriminator method from parent class. Instead returns
        CNN based discriminator to use with MNIST data.

        Returns:
            keras.Model: Critic
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
        """Train the critic on one batch of data

        Args:
            labels (tf.Tensor): One-hot encoded image labels
            images (tf.Tensor): Images

        Returns:
            tf.Tensor: Critic loss for the step
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
        """Train the generator on one batch of data

        Args:
            labels (tf.Tensor): One-hot encoded image labels

        Returns:
            tf.Tensor: Generator loss for the step
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
    """Subclass of the Trainer class, used to train a cWGAN on MNIST data"""

    def __init__(self, params_dict):
        """Constructor

        Args:
            params_dict (dict): Contains training/model parameters
        """
        self.save_dir = file_utils.make_save_directory("cWGAN")
        self.num_critic_iters = params_dict["num_critic_iters"]
        self.batch_size = params_dict["batch_size"]
        clip_value = params_dict["clip_value"]
        noise_dims = params_dict["noise_dims"]
        gen_lr = params_dict["gen_lr"]
        critic_lr = params_dict["critic_lr"]
        gp_weight = params_dict["gp_weight"]
        self.model = cWGAN_mnist(clip_value, noise_dims, gen_lr, critic_lr, gp_weight)

        self.data = data_utils.load_mnist_data()
        self.epochs = params_dict["epochs"]
        self.num_training_examples = len(self.data[0])
        self.weight_saving_interval = params_dict["weight_saving_interval"]

        self.critic_losses = []
        self.generator_losses = []
        self.wass_estimates = []

    def take_generator_step(self):
        """Override function from parent class in order to handel the image - label
        concatenation necessary for the critic.
        """

        labels, images = self.sample_batch_of_data()
        concat_real = data_utils.concatenate_images_labels(images, labels)
        generator_loss = self.model.train_generator(labels)
        predicted_images = self.model.make_generator_predictions(labels)
        concat_fake = data_utils.concatenate_images_labels(predicted_images, labels)
        self.generator_losses.append(generator_loss)
        real_output = self.model.critic(concat_real, training=False)
        fake_output = self.model.critic(concat_fake, training=False)
        wass_estimate = -self.model.critic_loss(real_output, fake_output)
        self.wass_estimates.append(wass_estimate)
