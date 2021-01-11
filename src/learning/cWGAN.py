import tensorflow as tf
from tensorflow import keras
import numpy as np
from tensorflow.keras.constraints import max_norm


class cWGAN:
    """
    Class implementing a conditional wasserstein generative adversarial
    network
    """

    def __init__(self, clip_value=0.01, noise_dims=4):
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
        pJet = keras.Input(shape=(4,), name="pjetIn")

        x = keras.layers.Dense(32, activation="relu", name="genx1")(pJet)
        x = keras.layers.Dense(64, activation="relu", name="genx2")(x)
        x = keras.layers.Dense(64, activation="relu", name="genx3")(x)

        y = keras.layers.Dense(32, activation="relu", name="geny1")(noise)
        y = keras.layers.Dense(64, activation="relu", name="geny2")(y)
        y = keras.layers.Dense(64, activation="relu", name="geny3")(y)

        concat = keras.layers.concatenate([x, y], name="concat")
        out = keras.layers.Dense(64, activation="relu", name="both1")(concat)
        out = keras.layers.Dense(64, activation="relu", name="both2")(concat)
        out = keras.layers.Dense(4, name="out")(out)
        return keras.Model([pJet, noise], out)

    def build_critic(self):
        pJet = keras.Input(shape=(4,))
        rJet = keras.Input(shape=(4,))

        x = keras.layers.Dense(32, activation="relu")(pJet)
        x = keras.layers.Dense(64, activation="relu")(x)
        x = keras.layers.Dense(64, activation="relu")(x)

        y = keras.layers.Dense(32, activation="relu")(rJet)
        y = keras.layers.Dense(64, activation="relu")(y)
        y = keras.layers.Dense(64, activation="relu")(y)

        concat = keras.layers.concatenate([x, y])
        out = keras.layers.Dense(128, activation="relu")(concat)
        out = keras.layers.Dense(128, activation="relu")(concat)
        out = keras.layers.Dense(1)(out)
        return keras.Model([pJet, rJet], out)

    @tf.function
    def critic_loss(self, real_output, fake_output):
        """
        The negative of the estimate of the wasserstein distance (negative
        because we want to perform gradient ascent on the critic)
        real_output - output of discriminator when given parton jets matched
        with real reco jets
        fake_output - output of discriminator when given parton jets and reco
        jets from the generator
        returns - wasserstein loss
        """
        loss = -(tf.math.reduce_mean(real_output) - tf.math.reduce_mean(fake_output))
        return loss

    @tf.function
    def generator_loss(self, fake_output):
        """
        Estimate of wasserstein loss for the generator
        fake_output - output of discriminator when given parton jets and reco
        jets from the generator
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
    def train_critic(self, pJets, rJets):
        """
        Train critic on one batch of data
        pJets - batch of parton jet 4-momenta
        rJets - batch of maching reco jet 4-momenta
        returns - the critic loss for the batch
        """
        noise = tf.random.uniform(
            (tf.shape(pJets)[0], self.noise_dims), 0, 1, tf.float32
        )
        with tf.GradientTape(persistent=True) as tape:
            generated_rJets = self.generator([pJets, noise], training=False)
            real_output = self.critic([pJets, rJets], training=True)
            fake_output = self.critic([pJets, generated_rJets], training=True)

            critic_loss_val = self.critic_loss(real_output, fake_output)

        critic_grads = tape.gradient(critic_loss_val, self.critic.trainable_variables)

        self.critic_optimizer.apply_gradients(
            zip(critic_grads, self.critic.trainable_variables)
        )

        return critic_loss_val

    @tf.function
    def train_generator(self, pJets):
        """
        Train generator on one batch of data
        pJets - batch of parton jet 4-momenta
        """

        noise = tf.random.uniform(
            (tf.shape(pJets)[0], self.noise_dims), 0, 1, tf.float32
        )

        with tf.GradientTape() as tape:
            generated_rJets = self.generator([pJets, noise], training=True)
            fake_output = self.critic([pJets, generated_rJets], training=False)
            generator_loss_val = self.generator_loss(fake_output)

        generator_grads = tape.gradient(
            generator_loss_val, self.generator.trainable_variables
        )
        self.generator_optimizer.apply_gradients(
            zip(generator_grads, self.generator.trainable_variables)
        )

        return generator_loss_val

    @tf.function
    def make_generator_predictions(self, pJets):
        noise = tf.random.uniform(
            (tf.shape(pJets)[0], self.noise_dims), 0, 1, tf.float32
        )
        predictions = self.generator([pJets, noise], training=False)
        return predictions


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

        x = keras.layers.Dense(10, activation="relu")(number_input)
        x = keras.layers.Dense(32, activation="relu")(x)

        y = keras.layers.Dense(self.noise_dims, activation="relu")(noise)
        y = keras.layers.Dense(self.noise_dims, activation="relu")(y)

        concat = keras.layers.concatenate([x, y])
        out = keras.layers.Dense(7 * 7 * 256, activation="relu")(concat)

        out = keras.layers.Reshape((7, 7, 256))(out)
        out = keras.layers.Conv2DTranspose(
            128, (5, 5), strides=(1, 1), padding="same", use_bias=False
        )(out)
        out = keras.layers.LeakyReLU()(out)
        out = keras.layers.Conv2DTranspose(
            64, (5, 5), strides=(2, 2), padding="same", use_bias=False
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

        number_input = keras.Input(shape=(10,))
        image = keras.Input(shape=(28, 28, 1))

        x = keras.layers.Dense(10, activation="relu")(number_input)
        x = keras.layers.Dense(100, activation="relu")(x)

        y = keras.layers.Conv2D(128, (5, 5), strides=(2, 2), padding="same")(image)
        y = keras.layers.LeakyReLU()(y)
        y = keras.layers.Conv2D(128, (5, 5), strides=(2, 2), padding="same")(y)
        y = keras.layers.LeakyReLU()(y)
        y = keras.layers.Conv2D(256, (5, 5), strides=(2, 2), padding="same")(y)
        y = keras.layers.LeakyReLU()(y)

        y = keras.layers.Flatten()(y)

        concat = keras.layers.concatenate([x, y])
        out = keras.layers.Dense(100)(concat)
        out = keras.layers.Dense(100)(out)
        out = keras.layers.Dense(50)(out)
        out = keras.layers.Dense(1)(out)

        return keras.Model([number_input, image], out)


def main():

    net = cWGAN_mnist(5, 64, 100)
    net.print_network()


if __name__ == "__main__":
    main()
