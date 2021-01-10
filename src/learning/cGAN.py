import tensorflow as tf
from tensorflow import keras


class cGAN:
    """
    Class implementing a conditional generative adversarial network
    """

    def __init__(self, noise_dims):
        self.noise_dims = noise_dims

        self.generator = self.make_generator()
        self.discriminator = self.make_discriminator()

    def make_generator(self):
        noise = keras.Input(shape=(self.noise_dims,), name="noiseIn")
        pJet = keras.Input(shape=(4,), name="pjetIn")

        x = keras.layers.Dense(32, activation="relu", name="genx1")(pJet)
        x = keras.layers.Dense(32, activation="relu", name="genx2")(x)
        x = keras.layers.Dense(32, activation="relu", name="genx3")(x)

        y = keras.layers.Dense(32, activation="relu", name="geny1")(noise)
        y = keras.layers.Dense(32, activation="relu", name="geny2")(y)
        y = keras.layers.Dense(32, activation="relu", name="geny3")(y)

        concat = keras.layers.concatenate([x, y], name="concat")
        out = keras.layers.Dense(32, activation="relu", name="both1")(concat)
        out = keras.layers.Dense(32, activation="relu", name="both2")(concat)
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
        """
        Calculate binary crossentropy loss for the descriminator
        real_output - output from discriminator after being given parton
                    jets matched with real reco jets
        fake_output - output from discriminator after being given parton
                    jets and reco jets from the generator
        returns - loss for the disciminator
        """
        real_loss = tf.keras.losses.BinaryCrossentropy(
            tf.ones_like(real_output), real_output
        )
        fake_loss = tf.keras.losses.BinaryCrossentropy(
            tf.zeros_like(fake_output), fake_output
        )
        total_loss = real_loss + fake_loss
        return total_loss

    def generator_loss(self, fake_output):
        """
        Calculate binary crossentropy loss for the generator
        fake_output - output from discriminator after being given parton
                    jets and reco jets from the generator
        returns - loss for the generator
        """
        return tf.keras.losses.BinaryCrossentropy(
            tf.ones_like(fake_output), fake_output
        )


def main():
    pass


if __name__ == "__main__":
    main()
